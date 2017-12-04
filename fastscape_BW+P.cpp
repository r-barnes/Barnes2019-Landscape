//This file contains a serial implementation of Braun and Willett's FastScape
//algorithm. The implementation was developed by adapting Fortran code provided
//by Braun and attempts to be a faithful reproduction of the ideas therein.
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
  fout<<"ncols "<<width<<"\n";
  fout<<"nrows "<<height<<"\n";
  fout<<"xllcorner 637500.000\n"; //Arbitrarily chosen value
  fout<<"yllcorner 206000.000\n"; //Arbitrarily chosen value
  fout<<"cellsize 500.000\n";     //Arbitrarily chosen value
  fout<<"NODATA_value -9999\n";
  for(int y=0;y<height;y++){
    for(int x=0;x<width;x++)
      fout<<h[y*width+x]<<" ";
    fout<<"\n";
  }
}



class FastScape_BW {
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
  int width;
  int height;
  int size;

  double *h;
  double *accum;
  int    *rec;
  int    *ndon;
  int    *stack;
  int    *donor;
  int    nshift[8];

  std::vector<int>    stack_start;

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
  FastScape_BW(const int width0, const int height0)
    : nshift{-1,-width0-1,-width0,-width0+1,1,width0+1,width0,width0-1}
  {
    Tmr_Overall.start();
    Tmr_Step1_Initialize.start();
    width  = width0;
    height = height0;
    size   = width*height;

    //Each edge node of the DEM is the root of a stack. The `stack_start` array
    //indicates where these nodes are located in the global stack. This is then
    //used to induce parallelism. Note that, especially at the beginning of the
    //program, there may be interior nodes with in depressions or with no local
    //gradient. These may become the starts of stacks as well, so the
    //`stack_start` array may be larger than just the perimeter of the DEM.
    stack_start.reserve(3*width+2*height);

    h      = new double[size];

    GenerateRandomTerrain();

    Tmr_Step1_Initialize.stop();
    Tmr_Overall.stop();
  }

  ~FastScape_BW(){
    delete[] h;
  }


 private:
  void ComputeReceivers(){
    //! computing receiver array
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
      const int c      = y*width+x;
      double max_slope = -DINFTY;
      int    max_n     = 0;
      for(int n=0;n<8;n++){
        double slope = (h[c] - h[c+nshift[n]])/dr[n];
        if(slope>max_slope){
          max_slope = slope;
          max_n     = n;
        }
      }
      if(max_slope>-DINFTY)
        rec[c] = max_n;
      else
        rec[c] = c;
    }   
  }


  void ComputeDonors(){
    //! initialising number of donors per node to 0
    for(int i=0;i<size;i++)
      ndon[i] = 0;

    //! computing donor arrays
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW)
        continue;
      const int n        = c+nshift[rec[c]];
      donor[8*n+ndon[n]] = c;
      ndon[n]++;
    }
  }


  void FindStack(const int c, int &nstack){
    for(int k=0;k<ndon[c];k++){
      int n           = donor[8*c+k];
      stack[nstack++] = n;
      FindStack(n,nstack);
    }
  }


  void GenerateStack(){
    //The `stack_start` array has an unpredictable size, so we fill it
    //dynamically and reset it each time. This should have minimal impact on the
    //algorithm's speed since std::vector's memory is not actually reallocated.    
    stack_start.clear();

    int nstack=0;
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW){
        stack_start.push_back(nstack);
        stack[nstack++] = c;
        FindStack(c,nstack);
      }
    }  
    //We add an additional note to the end of `stack_start` that serves as an
    //upper bound on the locations of the cells in the final stack. See b    
    stack_start.push_back(nstack);
  }


  void ComputeDraingeArea(){
    //! computing drainage area
    for(int i=0;i<size;i++)
      accum[i] = cell_area;

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
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
      const int c = y*width+x;
      h[c] += ueq*dt;
    }    
  }


  void StackErode(){
    #pragma omp parallel for schedule(dynamic)
    for(int ss=0;ss<stack_start.size()-1;ss++){
      //Execute the inner for loop as a task
      const int sstart = stack_start.at(ss);
      const int send   = stack_start.at(ss+1);
      for(int s=sstart;s<send;s++){
        const int c = stack[s];           //Cell from which flow originates
        if(rec[c]==NO_FLOW)
          continue;
        const int n = c+nshift[rec[c]];   //Cell receiving the flow

        const double fact = keq*dt*std::pow(accum[c],meq)/std::pow(dr[rec[c]],neq);
        const double hn   = h[n];
        const double h0   = h[c];
        double hnew       = h0;
        double hp         = h0;
        double diff       = 2*tol;
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

    //! initializing rec
    for(int i=0;i<size;i++)
      rec[i] = NO_FLOW;

    Tmr_Step1_Initialize.stop();

    for(int step=0;step<=nstep;step++){
      Tmr_Step2_DetermineReceivers.start ();   ComputeReceivers  (); Tmr_Step2_DetermineReceivers.stop ();
      Tmr_Step3_DetermineDonors.start    ();   ComputeDonors     (); Tmr_Step3_DetermineDonors.stop    ();
      Tmr_Step4_GenerateStack.start      ();   GenerateStack     (); Tmr_Step4_GenerateStack.stop      ();
      Tmr_Step5_FlowAcc.start            ();   ComputeDraingeArea(); Tmr_Step5_FlowAcc.stop            ();
      Tmr_Step6_Uplift.start             ();   AddUplift         (); Tmr_Step6_Uplift.stop             ();
      Tmr_Step7_Erosion.start            ();   StackErode        (); Tmr_Step7_Erosion.stop            ();

      if( step%20==0 )
        std::cout<<step<<std::endl;
    }

    delete[] accum;
    delete[] rec;
    delete[] ndon;
    delete[] stack;
    delete[] donor;

    Tmr_Overall.stop();

    std::cout<<"t Step1: Initialize         = "<<std::setw(15)<<Tmr_Step1_Initialize.elapsed()         <<" microseconds"<<std::endl;                 
    std::cout<<"t Step2: DetermineReceivers = "<<std::setw(15)<<Tmr_Step2_DetermineReceivers.elapsed() <<" microseconds"<<std::endl;                         
    std::cout<<"t Step3: DetermineDonors    = "<<std::setw(15)<<Tmr_Step3_DetermineDonors.elapsed()    <<" microseconds"<<std::endl;                      
    std::cout<<"t Step4: GenerateStack      = "<<std::setw(15)<<Tmr_Step4_GenerateStack.elapsed()      <<" microseconds"<<std::endl;                    
    std::cout<<"t Step5: FlowAcc            = "<<std::setw(15)<<Tmr_Step5_FlowAcc.elapsed()            <<" microseconds"<<std::endl;              
    std::cout<<"t Step6: Uplift             = "<<std::setw(15)<<Tmr_Step6_Uplift.elapsed()             <<" microseconds"<<std::endl;             
    std::cout<<"t Step7: Erosion            = "<<std::setw(15)<<Tmr_Step7_Erosion.elapsed()            <<" microseconds"<<std::endl;              
    std::cout<<"t Overall                   = "<<std::setw(15)<<Tmr_Overall.elapsed()                  <<" microseconds"<<std::endl;        
  }


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
  FastScape_BW tm(width,height);
  tm.run(nstep);
  std::cout<<"Calculation time = "<<tmr.elapsed()<<std::endl;

  PrintDEM("out_BW_class.dem", tm.getH(), width, height);

  return 0;
}
