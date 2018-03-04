#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>  //Used for OpenMP run-time functions
#include "random.hpp"
#include <vector>
#include "CumulativeTimer.hpp"

//Used to handle situations in which OpenMP is not available
//(This scenario has not been extensively tested)
#ifndef _OPENMP
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
  #define omp_get_max_threads() 1
#endif


///This is a quick-and-dirty, zero-dependency function for saving the outputs of
///the model in ArcGIS ASCII DEM format (aka Arc/Info ASCII Grid, AAIGrid).
///Production code for experimentation should probably use GeoTIFF or a similar
///format as it will have a smaller file size and, thus, save quicker.
void PrintDEM(
  const std::string filename, 
  const double *const h,
  const int width,
  const int height
){
  std::ofstream fout(filename.c_str());
  //Since the outer ring of the dataset is a halo used for simplifying
  //neighbour-finding logic, we do not save it to the output here.
  fout<<"ncols "<<(width- 2)<<"\n";
  fout<<"nrows "<<(height-2)<<"\n";
  fout<<"xllcorner 637500.000\n"; //Arbitrarily chosen value
  fout<<"yllcorner 206000.000\n"; //Arbitrarily chosen value
  fout<<"cellsize 500.000\n";     //Arbitrarily chosen value
  fout<<"NODATA_value -9999\n";   //Value which is guaranteed not to correspond to an actual data value
  for(int y=1;y<height-1;y++){
    for(int x=1;x<width-1;x++)
      fout<<h[y*width+x]<<" ";
    fout<<"\n";
  }
}



///The entire model is contained in a handy class, which makes it easy to set up
///and solve many such models.
class FastScape_RBPQ {
 private:
  //Value used to indicate that a cell had no downhill neighbour and, thus, does
  //not flow anywhere.
  const int    NO_FLOW = -1;
  const double SQRT2   = 1.414213562373095048801688724209698078569671875376948; //Yup, this is overkill.


 public:
  //NOTE: Having these constants specified in the class rather than globally
  //results in a significant speed loss. However, it is better to have them here
  //under the assumption that they'd be dynamic in a real implementation.
  const double keq       = 2e-6;   //Stream power equation constant (coefficient)
  const double neq       = 2;      //Stream power equation constant (slope modifier)
  const double meq       = 0.8;    //Stream power equation constant (area modifier)
  const double ueq       = 2e-3;   //Rate of uplift
  const double dt        = 1000.;  //Timestep interval
  const double dr[8]     = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2}; //Distance between adjacent cell centers on a rectangular grid arbitrarily scale to cell edge lengths of 1
  const double tol       = 1e-3;   //Tolerance for Newton-Rhapson convergence while solving implicit Euler
  const double cell_area = 40000;  //Area of a single cell


 private:
  int width;        //Width of DEM
  int height;       //Height of DEM
  int size;         //Size of DEM (width*height)

  //The dataset is set up so that the outermost edge is never actually used:
  //it's a halo which allows every cell that is actually processed to consider
  //its neighbours in all 8 directions without having to check first to see if
  //it is an edge cell. The ring of second-most outer cells is set to a fixed
  //value to which everything erodes (in this model).

  //Rec directions (also used for nshift offsets) - see below for details
  //1 2 3
  //0   4
  //7 6 5

  double *h;        //Digital elevation model (height)
  double *accum;    //Flow accumulation at each point
  int    *rec;      //Direction of receiving cell
  int    *donor;    //Indices of a cell's donor cells
  int    *ndon;     //How many donors a cell has
  int    *stack;    //Indices of cells in the order they should be processed
  int    nshift[8]; //Offset from a focal cell's index to its neighbours in terms of flat indexing

  //A level is a set of cells which can all be processed simultaneously.
  //Topologically, cells within a level are neither descendents or ancestors of
  //each other in a topological sorting, but are the same number of steps from
  //the edge of the dataset.
  int    *levels;   //Indices of locations in stack where a level begins and ends
  int    nlevel;    //Number of levels used

  int stack_width;  //Number of cells allowed in the stack
  int level_width;  //Number of cells allowed in a level

  CumulativeTimer Tmr_Step1_Initialize;
  CumulativeTimer Tmr_Step2_DetermineReceivers;
  CumulativeTimer Tmr_Step3_DetermineDonors;
  CumulativeTimer Tmr_Step4_GenerateOrder;
  CumulativeTimer Tmr_Step5_FlowAcc;
  CumulativeTimer Tmr_Step6_Uplift;
  CumulativeTimer Tmr_Step7_Erosion;
  CumulativeTimer Tmr_Overall;


 private:
  void GenerateRandomTerrain(){
    //srand(std::random_device()());
    for(int y=0;y<height;y++)
    for(int x=0;x<width;x++){
      const int c = y*width+x;
      h[c]  = rand()/(double)RAND_MAX;

      //Outer edge is set to 0 and never touched again. It is used only as a
      //convenience so we don't have to worry when a focal cell looks at its
      //neighbours.
      if(x == 0 || y==0 || x==width-1 || y==height-1)
        h[c] = 0;

      //Second outer-most edge is set to 0 and never touched again. This is the
      //baseline to which all cells would erode where it not for uplift. You can
      //think of this as being "sea level".
      if(x == 1 || y==1 || x==width-2 || y==height-2)
        h[c] = 0;
    }
  }  


 public:
  ///Initializing code
  FastScape_RBPQ(const int width0, const int height0)
    //Initialize code for finding neighbours of a cell
    : nshift{-1,-width0-1,-width0,-width0+1,1,width0+1,width0,width0-1}
  {
    Tmr_Overall.start();
    Tmr_Step1_Initialize.start();
    width  = width0;
    height = height0;
    size   = width*height;

    h      = new double[size];   //Memory for terrain height

    GenerateRandomTerrain();     //Could replace this with custom initializer

    Tmr_Step1_Initialize.stop();
    Tmr_Overall.stop();
  }



  ///Destructor: ensures that `h` is freed when the class goes out of scope
  ~FastScape_RBPQ(){
    delete[] h;
  }



 private:
  ///The receiver of a focal cell is the cell which receives the focal cells'
  ///flow. Here, we model the receiving cell as being the one connected to the
  ///focal cell by the steppest gradient. If there is no local gradient, than
  ///the special value NO_FLOW is assigned.
  void ComputeReceivers(){
    //! computing receiver array
    #pragma omp for simd collapse(2)
    for(int y=2;y<height-2;y++)
    for(int x=2;x<width-2;x++){
      const int c      = y*width+x;

      //The slope must be greater than zero for there to be downhill flow;
      //otherwise, the cell is marekd NO_FLOW
      double max_slope = 0;
      int    max_n     = NO_FLOW;

      for(int n=0;n<8;n++){
        double slope = (h[c] - h[c+nshift[n]])/dr[n];
        if(slope>max_slope){
          max_slope = slope;
          max_n     = n;
        }
      }
      rec[c] = max_n;
    }    
  }


  void ComputeDonors(){
    //The B&W method of developing the donor array has each focal cell F inform
    //its receiving cell R that F is a donor of R. Unfortunately, parallelizing
    //this is difficult because more than one cell might be informing R at any
    //given time. Atomics are a solution, but they impose a performance cost
    //(though using the latest and greatest hardware decreases this penalty).

    //Instead, we invert the operation. Each focal cell now examines its
    //neighbours to see if it receives from them. Each focal cell is then
    //guaranteed to have sole write-access to its location in the donor array.
    #pragma omp for simd collapse(2) schedule(static) nowait
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
      const int c = y*width+x;
      ndon[c] = 0;
      for(int ni=0;ni<8;ni++){
        const int n = c+nshift[ni];
        if(rec[n]!=NO_FLOW && n+nshift[rec[n]]==c){
          donor[8*c+ndon[c]] = n;
          ndon[c]++;
        }
      }
    }
  }

  void GenerateOrder(
    int *const stack,
    int *const levels,
    int &nlevel
  ){
    int nstack = 0;

    levels[0]  = 0;
    nlevel     = 1;

    //Outer edge
    #pragma omp for schedule(static) nowait
    for(int y=1;y<height-1;y++){
      stack[nstack++] = y*width+1;          assert(nstack<stack_width);
      stack[nstack++] = y*width+(width-2);  assert(nstack<stack_width);
    }

    #pragma omp for schedule(static) nowait
    for(int x=2;x<width-2;x++){
      stack[nstack++] =          1*width+x; assert(nstack<stack_width);
      stack[nstack++] = (height-2)*width+x; assert(nstack<stack_width);
    }

    //End of outer edge
    levels[nlevel++] = nstack; //Last cell of this level

    //Interior cells
    //TODO: Outside edge is always NO_FLOW. Maybe this can get loaded once?
    //Load cells without dependencies into the queue
    //TODO: Why can't I use nowait here?
    #pragma omp for collapse(2) schedule(static) 
    for(int y=2;y<height-2;y++)
    for(int x=2;x<width -2;x++){
      const int c = y*width+x;
      if(rec[c]==NO_FLOW){
        stack[nstack++] = c;                
        assert(nstack<stack_width);
      }
    }
    //Last cell of this level
    levels[nlevel++] = nstack;              
    assert(nlevel<level_width); 

    int level_bottom = -1;
    int level_top    = 0;

    #pragma omp barrier //Must ensure ComputeDonors is done before we start the following

    while(level_bottom<level_top){
      level_bottom = level_top;
      level_top    = nstack;
      for(int si=level_bottom;si<level_top;si++){
        const auto c = stack[si];
        for(int k=0;k<ndon[c];k++){
          const auto n    = donor[8*c+k];
          stack[nstack++] = n;              
          assert(nstack<stack_width); //TODO `<=` or `<` ?
        }
      }

      levels[nlevel++] = nstack; //Starting a new level      
    }

    //End condition for the loop places two identical entries
    //at the end of the stack. Remove one.
    nlevel--;
  }


  void ComputeFlowAcc(
    const int *const stack,
    const int *const levels,
    const int &nlevel
  ){
    for(int i=levels[0];i<levels[nlevel-1];i++){
      const int c = stack[i];
      accum[c] = cell_area;
    }

    for(int li=nlevel-3;li>=1;li--){
      const int lvlstart = levels[li];
      const int lvlend   = levels[li+1];
      for(int si=lvlstart;si<lvlend;si++){
        const int c = stack[si];
        for(int k=0;k<ndon[c];k++){
          const auto n = donor[8*c+k];
          accum[c]    += accum[n];
        }
      }
    }    
  }


  void AddUplift(
    const int *const stack,
    const int *const levels,
    const int &nlevel
  ){
    //Start at levels[1] so we don't elevate the outer edge
    #pragma omp simd
    for(int i=levels[1];i<levels[nlevel-1];i++){
      const int c = stack[i];
      h[c]       += ueq*dt; 
    }
  }


  void Erode(
    const int *const stack,
    const int *const levels,
    const int nlevel
  ){
    //#pragma omp parallel default(none)
    for(int li=2;li<nlevel-1;li++){
      const int lvlstart = levels[li];
      const int lvlend   = levels[li+1];
      #pragma omp simd
      for(int si=lvlstart;si<lvlend;si++){
        const int c = stack[si];          //Cell from which flow originates
        const int n = c+nshift[rec[c]];   //Cell receiving the flow

        const double length = dr[rec[c]];
        const double fact   = keq*dt*std::pow(accum[c],meq)/std::pow(length,neq);
        const double h0     = h[c];        //Elevation of focal cell
        const double hn     = h[n];        //Elevation of neighbouring (receiving, lower) cell
        double hnew         = h0;          //Current updated value of focal cell
        double hp           = h0;          //Previous updated value of focal cell
        double diff         = 2*tol;       //Difference between current and previous updated values
        while(std::abs(diff)>tol){
          hnew -= (hnew-h0+fact*std::pow(hnew-hn,neq))/(1.+fact*neq*std::pow(hnew-hn,neq-1));
          diff  = hnew - hp;
          hp    = hnew;
        }
        h[c] = hnew;
      }
    }
  }


 public:
  void run(const int nstep){
    Tmr_Overall.start();

    Tmr_Step1_Initialize.start();

    accum  = new double[size];
    rec    = new int[size];
    ndon   = new int[size];
    donor  = new int[8*size];

    //! initializing rec
    #pragma omp parallel for
    for(int i=0;i<size;i++)
      rec[i] = NO_FLOW;

    #pragma omp parallel for
    for(int i=0;i<size;i++)
      ndon[i] = 0;

    //TODO: Make smaller, explain max
    stack_width = std::max(300000,5*size/omp_get_max_threads()); //Number of stack entries available to each thread
    level_width = std::max(1000,size/omp_get_max_threads());   //Number of level entries available to each thread

    #pragma omp parallel
    {
      int *stack  = new int[stack_width]; //Indices of cells in the order they should be processed

      //It's difficult to know how much memory should be allocated for levels. For
      //a square DEM with isotropic dispersion this is approximately sqrt(E/2). A
      //diagonally tilted surface with isotropic dispersion may have sqrt(E)
      //levels. A tortorously sinuous river may have up to E*E levels. We
      //compromise and choose a number of levels equal to the perimiter because
      //why not?
      int *levels = new int[level_width]; //TODO
      int  nlevel = 0;

      Tmr_Step1_Initialize.stop();

      for(int step=0;step<=nstep;step++){
        Tmr_Step2_DetermineReceivers.start ();   ComputeReceivers  ();                      Tmr_Step2_DetermineReceivers.stop ();
        Tmr_Step3_DetermineDonors.start    ();   ComputeDonors     ();                      Tmr_Step3_DetermineDonors.stop    ();
        Tmr_Step4_GenerateOrder.start      ();   GenerateOrder     (stack,levels,nlevel);   Tmr_Step4_GenerateOrder.stop      ();
        Tmr_Step5_FlowAcc.start            ();   ComputeFlowAcc    (stack,levels,nlevel);   Tmr_Step5_FlowAcc.stop            ();
        Tmr_Step6_Uplift.start             ();   AddUplift         (stack,levels,nlevel);   Tmr_Step6_Uplift.stop             ();
        Tmr_Step7_Erosion.start            ();   Erode             (stack,levels,nlevel);   Tmr_Step7_Erosion.stop            ();
        #pragma omp barrier //Ensure threads synchronize after erosion so we calculate receivers correctly

        #pragma omp master
        if( step%20==0 )
          std::cout<<"p Step = "<<step<<std::endl;
      }

      delete[] stack;
      delete[] levels;
    }

    Tmr_Overall.stop();

    std::cout<<"t Step1: Initialize         = "<<std::setw(15)<<Tmr_Step1_Initialize.elapsed()         <<" microseconds"<<std::endl;                 
    std::cout<<"t Step2: DetermineReceivers = "<<std::setw(15)<<Tmr_Step2_DetermineReceivers.elapsed() <<" microseconds"<<std::endl;                         
    std::cout<<"t Step3: DetermineDonors    = "<<std::setw(15)<<Tmr_Step3_DetermineDonors.elapsed()    <<" microseconds"<<std::endl;                      
    std::cout<<"t Step4: GenerateOrder      = "<<std::setw(15)<<Tmr_Step4_GenerateOrder.elapsed()      <<" microseconds"<<std::endl;                    
    std::cout<<"t Step5: FlowAcc            = "<<std::setw(15)<<Tmr_Step5_FlowAcc.elapsed()            <<" microseconds"<<std::endl;              
    std::cout<<"t Step6: Uplift             = "<<std::setw(15)<<Tmr_Step6_Uplift.elapsed()             <<" microseconds"<<std::endl;             
    std::cout<<"t Step7: Erosion            = "<<std::setw(15)<<Tmr_Step7_Erosion.elapsed()            <<" microseconds"<<std::endl;              
    std::cout<<"t Overall                   = "<<std::setw(15)<<Tmr_Overall.elapsed()                  <<" microseconds"<<std::endl;        

    delete[] accum;
    delete[] rec;
    delete[] ndon;
    delete[] donor;
  }

      // std::cerr<<"Levels: ";
      // for(auto &l: levels)
      //   std::cerr<<l<<" ";
      // std::cerr<<std::endl;

  double* getH() const {
    return h;
  }
};







int main(int argc, char **argv){
  //feenableexcept(FE_ALL_EXCEPT);

  if(argc!=5){
    std::cerr<<"Syntax: "<<argv[0]<<" <Dimension> <Steps> <Output Name> <Seed>"<<std::endl;
    return -1;
  }

  seed_rand(std::stoul(argv[4]));

  std::cout<<"A FastScape RB+PQ"<<std::endl;
  std::cout<<"C Richard Barnes TODO"<<std::endl;
  std::cout<<"h git_hash    = "<<GIT_HASH<<std::endl;
  std::cout<<"m Random seed = "<<argv[4]<<std::endl;

  const int width  = std::stoi(argv[1]);
  const int height = std::stoi(argv[1]);
  const int nstep  = std::stoi(argv[2]);

  CumulativeTimer tmr(true);
  FastScape_RBPQ tm(width,height);
  tm.run(nstep);
  std::cout<<"t Total calculation time    = "<<std::setw(15)<<tmr.elapsed()<<" microseconds"<<std::endl;

  PrintDEM(argv[3], tm.getH(), width, height);

  return 0;
}
