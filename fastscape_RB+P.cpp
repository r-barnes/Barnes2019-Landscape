#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <iomanip>
#include "CumulativeTimer.hpp"
#include "Timer.hpp"



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



class FastScape_RBP {
 private:
  static constexpr double DINFTY  = std::numeric_limits<double>::infinity();

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
  int intsize;      //Size of DEM minus boundary `(width-2)*(height-2)`

  double *h;        //Digital elevation model (height)
  double *accum;    //Flow accumulation at each point
  int    *rec;      //Index of receiving cell
  int    *donor;    //Indices of a cell's donor cells
  int    *ndon;     //How many donors a cell has
  int    *stack;    //Indices of cells in the order they should be processed
  int    nshift[8]; //Offset from a focal cell's index to its neighbours

  //nshift offsets:
  //1 2 3
  //0   4
  //7 6 5

  int    *levels;   //Indices of locations in stack where a level begins and ends
  int    nlevel;    //Number of levels used
  
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
    }
  }  


 public:
  FastScape_RBP(const int width0, const int height0)
    : nshift{-1,-width0-1,-width0,-width0+1,1,width0+1,width0,width0-1}
  {
    Tmr_Overall.start();
    Tmr_Step1_Initialize.start();
    width  = width0;
    height = height0;
    size   = width*height;

    h      = new double[size];

    GenerateRandomTerrain();

    Tmr_Step1_Initialize.stop();
    Tmr_Overall.stop();
  }

  ~FastScape_RBP(){
    delete[] h;
  }


 private:
  void ComputeReceivers(){
    //! computing receiver array
    #pragma omp parallel for collapse(2)
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
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
    //! initialising number of donors per node to 0
    for(int i=0;i<size;i++)
      ndon[i] = 0;

    //! computing donor arrays
    // for(int c=0;c<size;c++){
    //   if(rec[c]==NO_FLOW)
    //     continue;
    //   const auto n       = c+nshift[rec[c]];
    //   donor[8*n+ndon[n]] = c;
    //   ndon[n]++;
    // }


    auto do_edge = [&](const int c){
      if(rec[c]==NO_FLOW)
        return;
      const auto n       = c+nshift[rec[c]];
      donor[8*n+ndon[n]] = c;
      ndon[n]++;
    };

    #pragma omp parallel for
    for(int y=1;y<height-1;y++){
      do_edge(y*width+0);
      do_edge(y*width+width-1);
    }

    #pragma omp parallel for
    for(int x=1;x<width-1;x++){
      do_edge(         0*width+x);
      do_edge((height-1)*width+x);
    }    

    //! computing donor arrays
    #pragma omp parallel for collapse(2)
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
      const int c = y*width+x;
      for(int ni=0;ni<8;ni++){
        const int n = c+nshift[ni];
        if(rec[n]==c){
          donor[8*c+ndon[n]] = c;
          ndon[n]++;
        }
      }
    }
  }

  void GenerateQueue(){
    int nstack = 0;
    int qpoint = 0;

    levels[0] = 0;
    nlevel    = 1;

    //Load cells without dependencies into the queue
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW)
        stack[nstack++] = c;
    }
    levels[nlevel++] = nstack; //Last cell of this level

    while(qpoint<nstack){
      const auto c = stack[qpoint];
      for(int k=0;k<ndon[c];k++){
        const auto n    = donor[8*c+k];
        stack[nstack++] = n;
      }

      //TODO: What's this about, then?
      qpoint++;
      if(qpoint==levels[nlevel-1] && nstack!=levels[nlevel-1])
        levels[nlevel++] = nstack; //Starting a new level      
    }
    // std::cerr<<"nstack final = "<<nstack<<std::endl;
  }

  void ComputeFlowAcc(){
    //! computing drainage area
    for(int i=0;i<size;i++)
      accum[i] = cell_area;

    // #pragma omp parallel default(none)
    // for(int li=nlevel-2;li>=0;li--){
    //   #pragma omp for
    //   for(int si=levels[li];si<levels[li+1];si++){
    //     const int c = stack[si];

    //     if(rec[c]!=NO_FLOW){
    //       const int n = c+nshift[rec[c]];
    //       accum[n]   += accum[c];
    //     }
    //   }
    // }    

    for(int s=size-1;s>=0;s--){
      const int c = stack[s];
      if(rec[c]!=NO_FLOW){
        const int n = c+nshift[rec[c]];
        accum[n]   += accum[c];
      }
    }    
  }


  void AddUplift(){
    //! adding uplift to landscape
    #pragma omp parallel for collapse(2)
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
      const int c = y*width+x;
      h[c] += ueq*dt;
    }
  }


  void Erode(){
    //#pragma omp parallel default(none)
    for(int li=0;li<nlevel-1;li++){
      const int lvlstart = levels[li];
      const int lvlend   = levels[li+1];
      const int lvlsize  = lvlend-lvlstart;
      #pragma omp parallel for if(lvlstart>2000)
      for(int si=levels[li];si<levels[li+1];si++){
        const int c = stack[si];          //Cell from which flow originates
        if(rec[c]==NO_FLOW)
          continue;
        const int n = c+nshift[rec[c]];   //Cell receiving the flow

        const double length = dr[rec[c]];
        const double fact   = keq*dt*std::pow(accum[c],meq)/std::pow(length,neq);
        const double h0     = h[c];      //Elevation of focal cell
        const double hn     = h[n];      //Elevation of neighbouring (receiving, lower) cell
        double hnew         = h0;        //Current updated value of focal cell
        double hp           = h0;        //Previous updated value of focal cell
        double diff         = 2*tol;     //Difference between current and previous updated values
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
    stack  = new int[size];
    donor  = new int[8*size];

    //It's difficult to know how much memory should be allocated for levels. For
    //a square DEM with isotropic dispersion this is approximately sqrt(E/2). A
    //diagonally tilted surface with isotropic dispersion may have sqrt(E)
    //levels. A tortorously sinuous river may have up to E*E levels. We
    //compromise and choose a number of levels equal to the perimiter because
    //why not?
    levels = new int[size]; //TODO: Make smaller to `2*width+2*height`

    //! initializing rec
    for(int i=0;i<size;i++)
      rec[i] = NO_FLOW;

    Tmr_Step1_Initialize.stop();

    for(int step=0;step<=nstep;step++){
      Tmr_Step2_DetermineReceivers.start ();   ComputeReceivers  (); Tmr_Step2_DetermineReceivers.stop ();
      Tmr_Step3_DetermineDonors.start    ();   ComputeDonors     (); Tmr_Step3_DetermineDonors.stop    ();
      Tmr_Step4_GenerateStack.start      ();   GenerateQueue     (); Tmr_Step4_GenerateStack.stop      ();
      Tmr_Step5_FlowAcc.start            ();   ComputeFlowAcc    (); Tmr_Step5_FlowAcc.stop            ();
      Tmr_Step6_Uplift.start             ();   AddUplift         (); Tmr_Step6_Uplift.stop             ();
      Tmr_Step7_Erosion.start            ();   Erode             (); Tmr_Step7_Erosion.stop            ();

      if( step%20==0 )
        std::cout<<step<<std::endl;
    }

    delete[] accum;
    delete[] rec;
    delete[] ndon;
    delete[] stack;
    delete[] donor;
    delete[] levels;

    Tmr_Overall.stop();

    std::cout<<"t Step1: Initialize         = "<<std::setw(15)<<Tmr_Step1_Initialize.elapsed()         <<" microseconds"<<std::endl;                 
    std::cout<<"t Step2: DetermineReceivers = "<<std::setw(15)<<Tmr_Step2_DetermineReceivers.elapsed() <<" microseconds"<<std::endl;                         
    std::cout<<"t Step3: DetermineDonors    = "<<std::setw(15)<<Tmr_Step3_DetermineDonors.elapsed()    <<" microseconds"<<std::endl;                      
    std::cout<<"t Step4: GenerateQueue      = "<<std::setw(15)<<Tmr_Step4_GenerateStack.elapsed()      <<" microseconds"<<std::endl;                    
    std::cout<<"t Step5: FlowAcc            = "<<std::setw(15)<<Tmr_Step5_FlowAcc.elapsed()            <<" microseconds"<<std::endl;              
    std::cout<<"t Step6: Uplift             = "<<std::setw(15)<<Tmr_Step6_Uplift.elapsed()             <<" microseconds"<<std::endl;             
    std::cout<<"t Step7: Erosion            = "<<std::setw(15)<<Tmr_Step7_Erosion.elapsed()            <<" microseconds"<<std::endl;              
    std::cout<<"t Overall                   = "<<std::setw(15)<<Tmr_Overall.elapsed()                  <<" microseconds"<<std::endl;        
  }

      // std::cerr<<"Levels: ";
      // for(auto &l: levels)
      //   std::cerr<<l<<" ";
      // std::cerr<<std::endl;

  double* getH() const {
    return h;
  }
};







int main(){
  //feenableexcept(FE_ALL_EXCEPT);

  const int width  = 501;
  const int height = 501;
  const int nstep  = 120;

  Timer tmr;
  FastScape_RBP tm(width,height);
  tm.run(nstep);
  std::cout<<"Calculation time = "<<tmr.elapsed()<<std::endl;

  PrintDEM("out_RB+P.dem", tm.getH(), width, height);

  return 0;
}
