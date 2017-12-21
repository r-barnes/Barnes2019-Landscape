#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <vector>
#include "CumulativeTimer.hpp"

#ifndef _OPENMP
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
  #define omp_get_max_threads() 1
#endif

void PrintDEM(
  const std::string filename, 
  const double *const h,
  const int width,
  const int height
){
  std::ofstream fout(filename.c_str());
  fout<<"ncols "<<(width- 2)<<"\n";
  fout<<"nrows "<<(height-2)<<"\n";
  fout<<"xllcorner 637500.000\n"; //Arbitrarily chosen value
  fout<<"yllcorner 206000.000\n"; //Arbitrarily chosen value
  fout<<"cellsize 500.000\n";     //Arbitrarily chosen value
  fout<<"NODATA_value -9999\n";
  for(int y=1;y<height-1;y++){
    for(int x=1;x<width-1;x++)
      fout<<h[y*width+x]<<" ";
    fout<<"\n";
  }
}



class FastScape_RBPQ {
 private:
  const int    NO_FLOW = -1;
  const double SQRT2   = 1.414213562373095048801688724209698078569671875376948;


 public:
  //NOTE: Having these constants specified in the class rather than globally
  //results in a significant speed loss. However, it is better to have them here
  //under the assumption that they'd be dynamic in a real implementation.
  const double keq       = 2e-6;
  const double neq       = 2;
  const double meq       = 0.8;
  const double ueq       = 2e-3;
  const double dt        = 1000.;
  const double dr[8]     = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};  
  const double tol       = 1e-3;
  const double cell_area = 40000;


 private:
  int width;        //Width of DEM
  int height;       //Height of DEM
  int size;         //Size of DEM (width*height)

  double *h;        //Digital elevation model (height)
  double *accum;    //Flow accumulation at each point
  int    *rec;      //Index of receiving cell
  int    *donor;    //Indices of a cell's donor cells
  int    *ndon;     //How many donors a cell has
  int    *stack;    //Indices of cells in the order they should be processed
  int    nshift[8]; //Offset from a focal cell's index to its neighbours

  int    t_stack_width; //Number of stack entries available to each thread
  int    t_level_width; //Number of level entries available to each thread

  //nshift offsets:
  //1 2 3
  //0   4
  //7 6 5

  int    *levels;   //Indices of locations in stack where a level begins and ends
  int    *nlevel;    //Number of levels used
  
  CumulativeTimer Tmr_Step1_Initialize;
  CumulativeTimer Tmr_Step2_DetermineReceivers;
  CumulativeTimer Tmr_Step3_DetermineDonors;
  CumulativeTimer Tmr_Step4_GenerateStack;
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
      if(x == 0 || y==0 || x==width-1 || y==height-1)
        h[c] = 0;
      if(x == 1 || y==1 || x==width-2 || y==height-2)
        h[c] = 0;
    }
  }  


 public:
  FastScape_RBPQ(const int width0, const int height0)
    : nshift{-1,-width0-1,-width0,-width0+1,1,width0+1,width0,width0-1}
  {
    Tmr_Overall.start();
    Tmr_Step1_Initialize.start();
    width  = width0;
    height = height0;
    size   = width*height;

    h      = new double[size];

    t_stack_width = std::max(100,2*size/omp_get_max_threads()); //TODO: Explain
    t_level_width = std::max(100,size/omp_get_max_threads());   //TODO: Explain, make smaller

    GenerateRandomTerrain();

    Tmr_Step1_Initialize.stop();
    Tmr_Overall.stop();
  }

  ~FastScape_RBPQ(){
    delete[] h;
  }

  void printDiagnostic(std::string msg){
    return;
    std::cerr<<"\n#################\n"<<msg<<std::endl;

    std::cerr<<"idx: "<<std::endl;
    for(int y=0;y<height;y++){
      for(int x=0;x<width;x++){
        std::cerr<<std::setw(6)<<std::setprecision(3)<<h[y*width+x];
        std::cerr<<"| ";
      }
      std::cerr<<"\n";
    }

    std::cerr<<"idx: "<<std::endl;
    for(int y=0;y<height;y++){
      for(int x=0;x<width;x++){
        std::cerr<<std::setw(6)<<(y*width+x);
        std::cerr<<"| ";
      }
      std::cerr<<"\n";
    }

    std::cerr<<"Rec: "<<std::endl;
    std::cerr<<"NO_FLOW = "<<NO_FLOW<<std::endl;
    for(int y=0;y<height;y++){
      for(int x=0;x<width;x++){
        std::cerr<<std::setw(6)<<rec[y*width+x];
        std::cerr<<"| ";
      }
      std::cerr<<"\n";
    }    

    std::cerr<<"Donor: "<<std::endl;
    for(int x=0;x<width;x++)
      std::cerr<<std::setw(24)<<x<<"|";
    std::cerr<<std::endl;
    for(int y=0;y<height;y++){
      for(int x=0;x<width;x++){
        const int c = y*width+x;
        for(int ni=0;ni<8;ni++)
          std::cerr<<std::setw(3)<<donor[8*c+ni];
        std::cerr<<"|";
      }
      std::cerr<<"\n";
    }    

    std::cerr<<"ndon: "<<std::endl;
    for(int y=0;y<height;y++){
      for(int x=0;x<width;x++){
        std::cerr<<std::setw(6)<<ndon[y*width+x];
        std::cerr<<"| ";
      }
      std::cerr<<"\n";
    }        
  }

 private:
  void ComputeReceivers(){
    //! computing receiver array
    #pragma omp barrier

    #pragma omp for collapse(2) schedule(static) nowait
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

    #pragma omp barrier

    #pragma omp for collapse(2) schedule(static) nowait
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

  void GenerateQueue(){
    const int tnum = omp_get_thread_num();
    auto tstack    = &stack[tnum*t_stack_width];
    auto tlevels   = &levels[tnum*t_level_width];
    int &tnlevel   = nlevel[tnum];

    //TODO: These loops are just a safety feature
    for(int i=0;i<t_stack_width;i++)
      tstack[i]  = -1;
    for(int i=0;i<t_level_width;i++)
      tlevels[i] = -1;

    int tnstack = 0; //Thread local
    tlevels[0] = 0;
    tnlevel    = 1;

    //Outer edge
    #pragma omp for schedule(static) nowait
    for(int y=1;y<height-1;y++){
      tstack[tnstack++] = y*width+1;          assert(tnstack<t_stack_width);
      tstack[tnstack++] = y*width+(width-2);  assert(tnstack<t_stack_width);
    }

    #pragma omp for schedule(static) nowait
    for(int x=1;x<width-1;x++){
      tstack[tnstack++] =          1*width+x; assert(tnstack<t_stack_width);
      tstack[tnstack++] = (height-2)*width+x; assert(tnstack<t_stack_width);
    }

    //End of outer edge
    tlevels[tnlevel++] = tnstack; //Last cell of this level

    //Interior cells
    //TODO: Outside edge is always NO_FLOW. Maybe this can get loaded once?
    //Load cells without dependencies into the queue
    //TODO: Why can't I use nowait here?
    #pragma omp for collapse(2) schedule(static) 
    for(int y=2;y<height-2;y++)
    for(int x=2;x<width -2;x++){
      const int c = y*width+x;
      if(rec[c]==NO_FLOW){
        tstack[tnstack++] = c;                assert(tnstack<t_stack_width);
      }
    }
    //Last cell of this level
    tlevels[tnlevel++] = tnstack;             assert(tnlevel<t_level_width); 

    int qpoint = 0;

    while(qpoint<tnstack){
      const auto c = tstack[qpoint];
      for(int k=0;k<ndon[c];k++){
        const auto n    = donor[8*c+k];
        tstack[tnstack++] = n;                assert(tnstack<t_stack_width);
      }

      //TODO: What's this about, then?
      qpoint++;
      if(qpoint==tlevels[tnlevel-1] && tnstack!=tlevels[tnlevel-1]){
        //Starting a new level      
        tlevels[tnlevel++] = tnstack;         assert(tnlevel<t_level_width); 
      }
    }
    // std::cerr<<"tnstack final = "<<tnstack<<std::endl;
  }


  void ComputeFlowAcc(){
    //! computing drainage area
    const int tnum = omp_get_thread_num();
    auto tstack    = &stack[tnum*t_stack_width];
    auto tlevels   = &levels[tnum*t_level_width];
    int &tnlevel   = nlevel[tnum];

    for(int i=tlevels[0];i<tlevels[tnlevel-1];i++){
      const int c = tstack[i];
      accum[c] = cell_area;
    }

    for(int li=tnlevel-2;li>=0;li--){
      for(int si=tlevels[li];si<tlevels[li+1];si++){
        const int c = tstack[si];

        if(rec[c]!=NO_FLOW){
          const int n = c+nshift[rec[c]];
          accum[n]   += accum[c];
        }
      }
    }    
  }


  void AddUplift(){
    const int tnum = omp_get_thread_num();
    auto tstack    = &stack[tnum*t_stack_width];
    auto tlevels   = &levels[tnum*t_level_width];
    int &tnlevel   = nlevel[tnum];

    //Start at tlevels[1] so we don't elevate the outer edge
    for(int i=tlevels[1];i<tlevels[tnlevel-1];i++){
      const int c = tstack[i];
      h[c]       += ueq*dt; 
    }
  }


  void Erode(){
    const int tnum = omp_get_thread_num();
    auto tstack    = &stack[tnum*t_stack_width];
    auto tlevels   = &levels[tnum*t_level_width];
    int &tnlevel   = nlevel[tnum];

    //#pragma omp parallel default(none)
    for(int li=0;li<tnlevel-1;li++){
      for(int si=tlevels[li];si<tlevels[li+1];si++){
        const int c = tstack[si];          //Cell from which flow originates
        if(rec[c]==NO_FLOW)
          continue;
        const int n = c+nshift[rec[c]];    //Cell receiving the flow

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
    stack  = new int[omp_get_max_threads()*t_stack_width]; //TODO: Explain
    nlevel = new int[omp_get_max_threads()];

    //It's difficult to know how much memory should be allocated for levels. For
    //a square DEM with isotropic dispersion this is approximately sqrt(E/2). A
    //diagonally tilted surface with isotropic dispersion may have sqrt(E)
    //levels. A tortorously sinuous river may have up to E*E levels. We
    //compromise and choose a number of levels equal to the perimiter because
    //why not?
    levels = new int[omp_get_max_threads()*t_level_width]; //TODO: Make smaller to `2*width+2*height`

    //! initializing rec
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0;i<size;i++)
        rec[i] = NO_FLOW;

      #pragma omp for
      for(int i=0;i<size;i++)
        ndon[i] = 0;
    }


    Tmr_Step1_Initialize.stop();

    #pragma omp parallel
    for(int step=0;step<=nstep;step++){
      Tmr_Step2_DetermineReceivers.start ();   ComputeReceivers  (); Tmr_Step2_DetermineReceivers.stop ();
      Tmr_Step3_DetermineDonors.start    ();   ComputeDonors     (); Tmr_Step3_DetermineDonors.stop    ();
      Tmr_Step4_GenerateStack.start      ();   GenerateQueue     (); Tmr_Step4_GenerateStack.stop      ();
      Tmr_Step5_FlowAcc.start            ();   ComputeFlowAcc    (); Tmr_Step5_FlowAcc.stop            ();
      Tmr_Step6_Uplift.start             ();   AddUplift         (); Tmr_Step6_Uplift.stop             ();
      Tmr_Step7_Erosion.start            ();   Erode             (); Tmr_Step7_Erosion.stop            ();

      #pragma omp master
      if( step%20==0 )
        std::cout<<"p Step = "<<step<<std::endl;
    }

    Tmr_Overall.stop();

    std::cout<<"t Step1: Initialize         = "<<std::setw(15)<<Tmr_Step1_Initialize.elapsed()         <<" microseconds"<<std::endl;                 
    std::cout<<"t Step2: DetermineReceivers = "<<std::setw(15)<<Tmr_Step2_DetermineReceivers.elapsed() <<" microseconds"<<std::endl;                         
    std::cout<<"t Step3: DetermineDonors    = "<<std::setw(15)<<Tmr_Step3_DetermineDonors.elapsed()    <<" microseconds"<<std::endl;                      
    std::cout<<"t Step4: GenerateOrder      = "<<std::setw(15)<<Tmr_Step4_GenerateStack.elapsed()      <<" microseconds"<<std::endl;                    
    std::cout<<"t Step5: FlowAcc            = "<<std::setw(15)<<Tmr_Step5_FlowAcc.elapsed()            <<" microseconds"<<std::endl;              
    std::cout<<"t Step6: Uplift             = "<<std::setw(15)<<Tmr_Step6_Uplift.elapsed()             <<" microseconds"<<std::endl;             
    std::cout<<"t Step7: Erosion            = "<<std::setw(15)<<Tmr_Step7_Erosion.elapsed()            <<" microseconds"<<std::endl;              
    std::cout<<"t Overall                   = "<<std::setw(15)<<Tmr_Overall.elapsed()                  <<" microseconds"<<std::endl;        

    delete[] accum;
    delete[] rec;
    delete[] ndon;
    delete[] stack;
    delete[] donor;
    delete[] levels;
    delete[] nlevel;    
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

  if(argc!=3){
    std::cerr<<"Syntax: "<<argv[0]<<" <Dimension> <Steps>"<<std::endl;
    return -1;
  }

  std::cout<<"A FastScape RB+PQ"<<std::endl;
  std::cout<<"C Richard Barnes TODO"<<std::endl;

  const int width  = std::stoi(argv[1]);
  const int height = std::stoi(argv[1]);
  const int nstep  = std::stoi(argv[2]);

  CumulativeTimer tmr(true);
  FastScape_RBPQ tm(width,height);
  tm.run(nstep);
  std::cout<<"t Total calculation time    = "<<std::setw(15)<<tmr.elapsed()<<" microseconds"<<std::endl;

  PrintDEM("out_RB+PQ.dem", tm.getH(), width, height);

  return 0;
}
