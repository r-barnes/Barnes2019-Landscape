#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include "random.hpp"
#include <vector>
#include "CumulativeTimer.hpp"



///This is a quick-and-dirty, zero-dependency function for saving the outputs of
///the model in ArcGIS ASCII DEM format (aka Arc/Info ASCII Grid, AAIGrid).
///Production code for experimentation should probably use GeoTIFF or a similar
///format as it will have a smaller file size and, thus, save quicker.
void PrintDEM(
  const std::string filename,
  const std::vector<double>& h,
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
class FastScape_RB {
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

  std::vector<double> h;        //Digital elevation model (height)
  std::vector<double> accum;    //Flow accumulation at each point
  std::vector<int>    rec;      //Direction of receiving cell
  std::vector<int>    donor;    //Indices of a cell's donor cells
  std::vector<int>    ndon;     //How many donors a cell has
  std::vector<int>    stack;    //Indices of cells in the order they should be processed
  std::array<int,8>   nshift;   //Offset from a focal cell's index to its neighbours in terms of flat indexing

  //A level is a set of cells which can all be processed simultaneously.
  //Topologically, cells within a level are neither descendents or ancestors of
  //each other in a topological sorting, but are the same number of steps from
  //the edge of the dataset.
  std::vector<int>    levels;   //Indices of locations in stack where a level begins and ends
  int    nlevel;    //Number of levels used

  //Timers for keeping track of how long each part of the code takes
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
      h[c]  = uniform_rand_real(0,1);

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
  FastScape_RB(const int width0, const int height0)
    //Initialize code for finding neighbours of a cell
    : nshift{-1,-width0-1,-width0,-width0+1,1,width0+1,width0,width0-1}
  {
    Tmr_Overall.start();
    Tmr_Step1_Initialize.start();
    width  = width0;
    height = height0;
    size   = width*height;

    h.resize(size);   //Memory for terrain height

    GenerateRandomTerrain();     //Could replace this with custom initializer

    Tmr_Step1_Initialize.stop();
    Tmr_Overall.stop();
  }



 private:
  ///The receiver of a focal cell is the cell which receives the focal cells'
  ///flow. Here, we model the receiving cell as being the one connected to the
  ///focal cell by the steppest gradient. If there is no local gradient, than
  ///the special value NO_FLOW is assigned.
  void ComputeReceivers(){
    //Edge cells do not have receivers because they do not distribute their flow
    //to anywhere.
    for(int y=2;y<height-2;y++)
    for(int x=2;x<width-2;x++){
      const int c      = y*width+x;

      //The slope must be greater than zero for there to be downhill flow;
      //otherwise, the cell is marked NO_FLOW.
      double max_slope = 0;        //Maximum slope seen so far amongst neighbours
      int    max_n     = NO_FLOW;  //Direction of neighbour which had maximum slope to focal cell

      //Loop over neighbours
      for(int n=0;n<8;n++){
        const double slope = (h[c] - h[c+nshift[n]])/dr[n]; //Slope to neighbour n
        if(slope>max_slope){    //Is this the steepest slope we've seen?
          max_slope = slope;    //If so, make a note of the slope
          max_n     = n;        //And which cell it came from
        }
      }
      rec[c] = max_n;           //Having considered all neighbours, this is the steepest
    }
  }



  ///The donors of a focal cell are the neighbours from which it receives flow.
  ///Here, we identify those neighbours by inverting the Receivers array.
  void ComputeDonors(){
    //Initially, we claim that each cell has no donors.
    for(int i=0;i<size;i++)
      ndon[i] = 0;

    //Looping across all cells
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW)
        continue;
      //If this cell passes flow to a downhill cell, make a note of it in that
      //downhill cell's donor array and increment its donor counter
      const auto n       = c+nshift[rec[c]];
      donor[8*n+ndon[n]] = c;
      ndon[n]++;
    }
  }



  ///Cells must be ordered so that they can be traversed such that higher cells
  ///are processed before their lower neighbouring cells. This method creates
  ///such an order. It also produces a list of "levels": cells which are,
  ///topologically, neither higher nor lower than each other. Cells in the same
  ///level can all be processed simultaneously without having to worry about
  ///race conditions.
  void GenerateOrder(){
    int nstack = 0;    //Number of cells currently in the stack

    //Since each value of the `levels` array is later used as the starting value
    //of a for-loop, we include a zero at the beginning of the array.
    levels[0] = 0;
    nlevel    = 1;     //Note that array now contains a single value

    //Load cells without dependencies into the queue. This will include all of
    //the edge cells.
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW)
        stack[nstack++] = c;
    }
    levels[nlevel++] = nstack; //Last cell of this level

    //Start with level_bottom=-1 so we get into the loop, it is immediately
    //replaced by level_top.
    int level_bottom = -1;         //First cell of the current level
    int level_top    =  0;         //Last cell of the current level

    while(level_bottom<level_top){ //Enusre we parse all the cells
      level_bottom = level_top;    //The top of the previous level we considered is the bottom of the current level
      level_top    = nstack;       //The new top is the end of the stack (last cell added from the previous level)
      for(int si=level_bottom;si<level_top;si++){
        const auto c = stack[si];
        //Load donating neighbours of focal cell into the stack
        for(int k=0;k<ndon[c];k++){
          const auto n = donor[8*c+k];
          stack[nstack++] = n;
        }
      }

      levels[nlevel++] = nstack; //Start a new level
    }

    //End condition for the loop places two identical entries
    //at the end of the stack. Remove one.
    nlevel--;

    assert(levels[nlevel-1]==nstack);
  }



  ///Compute the flow accumulation for each cell: the number of cells whose flow
  ///ultimately passes through the focal cell multiplied by the area of each
  ///cell. Each cell could also have its own weighting based on, say, average
  ///rainfall.
  void ComputeFlowAcc(){
    //Initialize cell areas to their weights. Here, all the weights are the
    //same.
    for(int i=0;i<size;i++)
      accum[i] = cell_area;

    //Highly-elevated cells pass their flow to less elevated neighbour cells.
    //The queue is ordered so that higher cells are keyed to higher indices in
    //the queue; therefore, parsing the queue in reverse ensures that fluid
    //flows downhill.
    for(int s=size-1;s>=0;s--){
      const int c = stack[s];
      if(rec[c]!=NO_FLOW){
        const int n = c+nshift[rec[c]];
        accum[n]   += accum[c];
      }
    }
  }



  ///Raise each cell in the landscape by some amount, otherwise it wil get worn
  ///flat (in this model, with these settings)
  void AddUplift(){
    //We exclude two exterior rings of cells in this example. The outermost ring
    //(the edges of the dataset) allows us to ignore the edges of the dataset,
    //the second-most outer ring (the cells bordering the edge cells of the
    //dataset) are fixed to a specified height in this model. All other cells
    //have heights which actively change and they are altered here.
    for(int y=2;y<height-2;y++)
    for(int x=2;x<width-2;x++){
      const int c = y*width+x;
      h[c] += ueq*dt;
    }
  }



  ///Decrease he height of cells according to the stream power equation; that
  ///is, based on a constant K, flow accumulation A, the local slope between
  ///the cell and its receiving neighbour, and some judiciously-chosen constants
  ///m and n.
  ///    h_next = h_current - K*dt*(A^m)*(Slope)^n
  ///We solve this equation implicitly to preserve accuracy
  void Erode(){
    for(int s=0;s<size;s++){
        const int c = stack[s];            //Cell from which flow originates
        if(rec[c]==NO_FLOW)              //Ignore cells with no receiving neighbour
          continue;
        const int n = c+nshift[rec[c]];  //Cell receiving the flow

        const double length = dr[rec[c]];
        //`fact` contains a set of values which are constant throughout the integration
        const double fact   = keq*dt*std::pow(accum[c],meq)/std::pow(length,neq);
        const double h0     = h[c];      //Elevation of focal cell
        const double hn     = h[n];      //Elevation of neighbouring (receiving, lower) cell
        double hnew         = h0;        //Current updated value of focal cell
        double hp           = h0;        //Previous updated value of focal cell
        double diff         = 2*tol;     //Difference between current and previous updated values
        while(std::abs(diff)>tol){       //Newton-Rhapson method (run until subsequent values differ by less than a tolerance, which can be set to any desired precision)
          hnew -= (hnew-h0+fact*std::pow(hnew-hn,neq))/(1.+fact*neq*std::pow(hnew-hn,neq-1));
          diff  = hnew - hp;             //Difference between previous and current value of the iteration
          hp    = hnew;                  //Update previous value to new value
        }
        h[c] = hnew;                     //Update value in array
    }
  }



 public:

  ///Run the model forward for a specified number of timesteps. No new
  ///initialization is done. This allows the model to be stopped, the terrain
  ///altered, and the model continued. For space-efficiency, a number of
  ///temporary arrays are created each time this is run, so repeatedly running
  ///this function for the same model will likely not be performant due to
  ///reallocations. If that is your use case, you'll want to modify your code
  ///appropriately.
  void run(const int nstep){
    Tmr_Overall.start();

    Tmr_Step1_Initialize.start();

    accum.resize(  size);  //Stores flow accumulation
    rec.resize  (  size);  //Array of Receiver directions
    ndon.resize (  size);  //Number of donors each cell has
    donor.resize(8*size);  //Array listing the donors of each cell (up to 8 for a rectangular grid)
    stack.resize(  size);  //Order in which to process cells

    //It's difficult to know how much memory should be allocated for levels. For
    //a square DEM with isotropic dispersion this is approximately sqrt(E/2). A
    //diagonally tilted surface with isotropic dispersion may have sqrt(E)
    //levels. A tortorously sinuous river may have up to E*E levels. We
    //compromise and choose a number of levels equal to the perimiter because
    //why not?
    levels.resize(2*width+2*height);

    ///All receivers initially point to nowhere
    for(int i=0;i<size;i++)
      rec[i] = NO_FLOW;

    Tmr_Step1_Initialize.stop();

    for(int step=0;step<=nstep;step++){
      Tmr_Step2_DetermineReceivers.start ();   ComputeReceivers  (); Tmr_Step2_DetermineReceivers.stop ();
      Tmr_Step3_DetermineDonors.start    ();   ComputeDonors     (); Tmr_Step3_DetermineDonors.stop    ();
      Tmr_Step4_GenerateOrder.start      ();   GenerateOrder     (); Tmr_Step4_GenerateOrder.stop      ();
      Tmr_Step5_FlowAcc.start            ();   ComputeFlowAcc    (); Tmr_Step5_FlowAcc.stop            ();
      Tmr_Step6_Uplift.start             ();   AddUplift         (); Tmr_Step6_Uplift.stop             ();
      Tmr_Step7_Erosion.start            ();   Erode             (); Tmr_Step7_Erosion.stop            ();

      if( step%20==0 ) //Show progress
        std::cout<<"p Step = "<<step<<std::endl;
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

    //Free up memory, except for the resulting landscape height field prior to
    //exiting so that unnecessary space is not used when the model is not being
    //run.
    accum .clear();   accum .shrink_to_fit();
    rec   .clear();   rec   .shrink_to_fit();
    ndon  .clear();   ndon  .shrink_to_fit();
    stack .clear();   stack .shrink_to_fit();
    donor .clear();   donor .shrink_to_fit();
    levels.clear();   levels.shrink_to_fit();
  }



  ///Returns a pointer to the data so that it can be copied, printed, &c.
  std::vector<double>& getH() {
    return h;
  }
};







int main(int argc, char **argv){
  //Enable this to stop the program if a floating-point exception happens
  //feenableexcept(FE_ALL_EXCEPT);

  if(argc!=5){
    std::cerr<<"Syntax: "<<argv[0]<<" <Dimension> <Steps> <Output Name> <Seed>"<<std::endl;
    return -1;
  }

  const int         width       = std::stoi (argv[1]);
  const int         height      = std::stoi (argv[1]);
  const int         nstep       = std::stoi (argv[2]);
  const std::string output_name =            argv[3] ;
  const auto        rand_seed   = std::stoul(argv[4]);

  seed_rand(rand_seed);

  //Uses the RichDEM machine-readable line prefixes
  //Name of algorithm
  std::cout<<"A FastScape RB"<<std::endl;
  //Citation for algorithm
  std::cout<<"C Richard Barnes TODO"<<std::endl;
  //Git hash of code used to produce outputs of algorithm
  std::cout<<"h git_hash    = "<<GIT_HASH<<std::endl;
  //Random seed used to produce outputs
  std::cout<<"m Random seed = "<<rand_seed<<std::endl;

  CumulativeTimer tmr(true);
  FastScape_RB tm(width,height);
  tm.run(nstep);
  std::cout<<"t Total calculation time    = "<<std::setw(15)<<tmr.elapsed()<<" microseconds"<<std::endl;

  PrintDEM(output_name, tm.getH(), width, height);

  return 0;
}
