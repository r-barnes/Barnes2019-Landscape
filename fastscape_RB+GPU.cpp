#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include "random.hpp"
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



class FastScape_RBGPU {
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

  //nshift offsets:
  //1 2 3
  //0   4
  //7 6 5
  int    nshift[8]; //Offset from a focal cell's index to its neighbours

  int    *levels;   //Indices of locations in stack where a level begins and ends
  int    nlevel;    //Number of levels used
  
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
      if(x == 0 || y==0 || x==width-1 || y==height-1)
        h[c] = 0;
      if(x == 1 || y==1 || x==width-2 || y==height-2)
        h[c] = 0;
    }
  }  


 public:
  FastScape_RBGPU(const int width0, const int height0)
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

  ~FastScape_RBGPU(){
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
    #pragma acc parallel loop independent gang collapse(2) present(this,h,rec,width,height,dr,nshift)
    for(int y=2;y<height-2;y++){
      #pragma acc loop independent
    for(int x=2;x<width-2;x++){
      const int c      = y*width+x;

      //The slope must be greater than zero for there to be downhill flow;
      //otherwise, the cell is marekd NO_FLOW
      double max_slope = 0;
      int    max_n     = NO_FLOW;

      #pragma acc loop seq
      for(int n=0;n<8;n++){
        double slope = (h[c] - h[c+nshift[n]])/dr[n];
        if(slope>max_slope){
          max_slope = slope;
          max_n     = n;
        }
      }
      rec[c] = max_n;
    }    }
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

    #pragma acc parallel loop collapse(2) present(this,rec,width,height,nshift,donor,ndon)
    for(int y=1;y<height-1;y++)
    #pragma acc loop independent
    for(int x=1;x<width-1;x++){
      const int c = y*width+x;
      ndon[c] = 0;
      #pragma acc loop seq
      for(int ni=0;ni<8;ni++){
        const int n = c+nshift[ni];
        if(rec[n]!=NO_FLOW && n+nshift[rec[n]]==c){
          donor[8*c+ndon[c]] = n;
          ndon[c]++;
        }
      }
    }
  }

  void GenerateOrder(){
    int nstack = 0;    //Number of levels used in the stack

    levels[0] = 0;
    nlevel    = 1;

    //TODO: Outside edge is always NO_FLOW. Maybe this can get loaded once?
    //Load cells without dependencies into the queue
    #pragma acc parallel loop present(this,rec,stack)
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW){
        int mystack;
        #pragma acc atomic capture
        mystack = nstack++;
        stack[mystack] = c;
      }
    }
    levels[nlevel++] = nstack; //Last cell of this level

    int level_bottom = 0;
    int level_top    = 0;

    //#pragma acc loop seq
    while(nstack<size){
      std::cerr<<nstack<<" - "<<size<<std::endl;
      level_bottom = level_top;
      level_top    = nstack;
      #pragma acc parallel loop independent copy(nstack) present(this,stack,donor,levels)
      for(int qpoint=level_bottom;qpoint<level_top;qpoint++){
        const auto c = stack[qpoint];
        #pragma acc loop seq
        for(int k=0;k<ndon[c];k++){
          const auto n    = donor[8*c+k];
          int mystack;
          #pragma acc atomic capture
          mystack = nstack++;
          stack[mystack] = n;
        }
      }

      //TODO: What's this about, then?
      if(nstack!=levels[nlevel-1])
        levels[nlevel++] = nstack; //Starting a new level      
    }
  }


  void ComputeFlowAcc(){
    //! computing drainage area
    #pragma acc parallel loop present(this,accum)
    for(int i=0;i<size;i++)
      accum[i] = cell_area;

    #pragma acc parallel loop seq present(this,accum,levels,nshift,rec)
    for(int li=nlevel-2;li>=0;li--){
      #pragma acc loop 
      for(int si=levels[li];si<levels[li+1];si++){
        const int c = stack[si];
        if(rec[c]!=NO_FLOW){
          const int n = c+nshift[rec[c]];
          accum[n]   += accum[c];
        }
      }
    }    
  }


  void AddUplift(){
    //! adding uplift to landscape
    #pragma acc parallel loop collapse(2) present(this,h)
    for(int y=2;y<height-2;y++)
    for(int x=2;x<width-2;x++){
      const int c = y*width+x;
      h[c]       += ueq*dt; 
    }
  }


  void Erode(){
    #pragma acc parallel loop seq present(this,levels,stack,nshift,rec,accum,h)
    for(int li=0;li<nlevel-1;li++){
      const int lvlstart = levels[li];
      const int lvlend   = levels[li+1];
      #pragma acc loop 
      for(int si=lvlstart;si<lvlend;si++){
        const int c = stack[si];          //Cell from which flow originates
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
        #pragma acc loop seq
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
    stack  = new int[size];

    //It's difficult to know how much memory should be allocated for levels. For
    //a square DEM with isotropic dispersion this is approximately sqrt(E/2). A
    //diagonally tilted surface with isotropic dispersion may have sqrt(E)
    //levels. A tortorously sinuous river may have up to E*E levels. We
    //compromise and choose a number of levels equal to the perimiter because
    //why not?
    levels = new int[size]; //TODO: Make smaller to `2*width+2*height`

    #pragma acc enter data copyin(this[0:1],h[0:size],nshift[0:8]) create(accum[0:size],rec[0:size],ndon[0:size],donor[0:8*size],stack[0:size],levels[0:size])

    //! initializing rec
    #pragma acc parallel loop present(rec)
    for(int i=0;i<size;i++)
      rec[i] = NO_FLOW;

    #pragma acc parallel loop present(ndon)
    for(int i=0;i<size;i++)
      ndon[i] = 0;


    Tmr_Step1_Initialize.stop();

    //#pragma acc kernels present(accum,rec,ndon,donor,stack,nlevel,levels)
    //#pragma acc loop seq
    for(int step=0;step<=nstep;step++){
      ComputeReceivers  ();
      ComputeDonors     ();
      GenerateOrder     ();
      ComputeFlowAcc    ();
      AddUplift         ();
      Erode             ();
    }

    Tmr_Overall.stop();

    std::cout<<"t Step1: Initialize         = "<<std::setw(15)<<Tmr_Step1_Initialize.elapsed()         <<" microseconds"<<std::endl;                 
    std::cout<<"t Overall                   = "<<std::setw(15)<<Tmr_Overall.elapsed()                  <<" microseconds"<<std::endl;        

    #pragma acc exit data copyout(h[0:size]) delete(this,accum,rec,ndon,donor,stack,nlevel,levels)

    delete[] accum;
    delete[] rec;
    delete[] ndon;
    delete[] stack;
    delete[] donor;
    delete[] levels;
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

  seed_rand(std::stoi(argv[4]));

  std::cout<<"A FastScape RB+GPU"<<std::endl;
  std::cout<<"C Richard Barnes TODO"<<std::endl;
  std::cout<<"h git_hash    = "<<GIT_HASH<<std::endl;
  std::cout<<"m Random seed = "<<argv[4]<<std::endl;

  const int width  = std::stoi(argv[1]);
  const int height = std::stoi(argv[1]);
  const int nstep  = std::stoi(argv[2]);

  CumulativeTimer tmr(true);
  FastScape_RBGPU tm(width,height);
  tm.run(nstep);
  std::cout<<"t Total calculation time    = "<<std::setw(15)<<tmr.elapsed()<<" microseconds"<<std::endl;

  PrintDEM(argv[3], tm.getH(), width, height);

  return 0;
}
