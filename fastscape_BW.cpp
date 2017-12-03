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

const double keq = 2e-6;
const double neq = 2;
const double meq = 0.8;
const double ueq = 2e-3;
const double dt  = 1000.;

constexpr double DINFTY  = std::numeric_limits<double>::infinity();
const     int    NO_FLOW = -1;

void find_stack(
  const int c, 
  const std::vector<int> &donor,
  const std::vector<int> &ndon,
  int SIZE, 
  std::vector<int> &stack,
  int &nstack
){
  for(int k=0;k<ndon[c];k++){
    int n           = donor[8*c+k];
    stack[nstack++] = n;
    find_stack(n,donor,ndon,SIZE,stack,nstack);
  }
}


void PrintDEM(
  const std::string filename, 
  const std::vector<double> &h,
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



int main(){
  CumulativeTimer Tmr_Step1_Initialize;
  CumulativeTimer Tmr_Step2_DetermineReceivers;
  CumulativeTimer Tmr_Step3_DetermineDonors;
  CumulativeTimer Tmr_Step4_GenerateStack;
  CumulativeTimer Tmr_Step5_FlowAcc;
  CumulativeTimer Tmr_Step6_Uplift;
  CumulativeTimer Tmr_Step7_Erosion;
  CumulativeTimer Tmr_Overall;

  Tmr_Overall.start();

  //feenableexcept(FE_ALL_EXCEPT);

  //! defining size of the problem
  const int WIDTH  = 501;
  const int HEIGHT = 501;
  const int SIZE   = WIDTH*HEIGHT;

  //!    allocating memory
  std::vector<double> h    (  SIZE);
  std::vector<double> accum(  SIZE);
  std::vector<double> length( SIZE);
  std::vector<int>    rec  (  SIZE);
  std::vector<int>    ndon (  SIZE);
  std::vector<int>    stack(  SIZE);
  std::vector<int>    donor(8*SIZE);
  

  //Neighbours
  const std::vector<int> nshift= {{-1,-WIDTH-1,-WIDTH,-WIDTH+1,1,WIDTH+1,WIDTH,WIDTH-1}};
  const double SQRT2  = 1.414213562373095048801688724209698078569671875376948;
  const double dr[8]  = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};

  //defining geometrical and temporal constants
  const double xl    = 100.e3;
  const double yl    = 100.e3;
  const double dx    = xl/(WIDTH-1);
  const double dy    = yl/(HEIGHT-1);
  const int    nstep = 120;
  const double tol   = 1.e-3;

  //! generating initial topography
  Tmr_Step1_Initialize.start();
  for(int y=0;y<HEIGHT;y++)
  for(int x=0;x<WIDTH;x++){
    const int c = y*WIDTH+x;
    h[c]  = rand()/(double)RAND_MAX;
    if(x == 0 || y==0 || x==WIDTH-1 || y==HEIGHT-1)
      h[c] = 0;
  }
  Tmr_Step1_Initialize.stop();  


  //! begining of time stepping
  for(int istep=0;istep<nstep;istep++){

    //! initializing rec and length
    Tmr_Step2_DetermineReceivers.start();
    for(int i=0;i<SIZE;i++){
      rec[i]    = i;
      length[i] = 0;
    }

    //! computing receiver array
    for(int y=1;y<HEIGHT-1;y++)
    for(int x=1;x<WIDTH-1;x++){
      const int c      = y*WIDTH+x;
      double max_slope = -DINFTY;
      int    max_n     = 0;
      double slope_n   = 0;
      for(int n=0;n<8;n++){
        double slope = (h[c] - h[c+nshift[n]])/dr[n];
        if(slope>max_slope){
          max_slope = slope;
          max_n     = n;
          slope_n   = dr[n];
        }
      }
      if(max_slope>-DINFTY){
        rec[c]    = c+nshift[max_n];
        length[c] = slope_n;
      }
    }
    Tmr_Step2_DetermineReceivers.stop();


    //! initialising number of donors per node to 0
    Tmr_Step3_DetermineDonors.start();    
    for(int i=0;i<SIZE;i++)
      ndon[i] = 0;

    //! computing donor arrays
    for(int c=0;c<SIZE;c++){
      if(rec[c]==c)
        continue;
      const int n        = rec[c];
      donor[8*n+ndon[n]] = c;
      ndon[n]++;
    }
    Tmr_Step3_DetermineDonors.stop();


    //! computing stack
    Tmr_Step4_GenerateStack.start();    
    int nstack=0;
    for(int c=0;c<SIZE;c++){
      if(rec[c]==c){
        stack[nstack++] = c;
        find_stack(c,donor,ndon,SIZE,stack,nstack);
      }
    }
    Tmr_Step4_GenerateStack.stop();


    //! computing drainage area
    Tmr_Step5_FlowAcc.start();    
    for(int i=0;i<SIZE;i++)
      accum[i] = dx*dy;

    for(int s=SIZE-1;s>=0;s--){
      const int c = stack[s];
      if(rec[c]!=c){
        const int n = rec[c];
        accum[n] += accum[c];
      }
    }
    Tmr_Step5_FlowAcc.stop();


    //! adding uplift to landscape
    Tmr_Step6_Uplift.start();    
    for(int y=1;y<HEIGHT-1;y++)
    for(int x=1;x<WIDTH-1;x++){
      int c = y*WIDTH+x;
      h[c] += ueq*dt;
    }
    Tmr_Step6_Uplift.stop();


    // std::cout<<"rec: ";
    // for(int c=3*WIDTH;c<4*WIDTH;c++){
    //   if(rec[c]==c)
    //     std::cout<<-1<<" ";
    //   else
    //     std::cout<<rec[c]<<" ";
    // }
    // std::cout<<std::endl;

    // std::cerr<<"length: ";
    // for(int i=3*WIDTH;i<4*WIDTH;i++)
    //   std::cerr<<length[i]<<" ";
    // std::cerr<<std::endl;

    //! computing erosion
    Tmr_Step7_Erosion.start();    
    for(int s=0;s<SIZE;s++){
      const int c = stack[s]; //Cell from which flow originates
      const int n = rec[c];   //Cell receiving the flow
      if(n==c)
        continue;
      const double fact = keq*dt*std::pow(accum[c],meq)/std::pow(length[c],neq);
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
    Tmr_Step7_Erosion.stop();

    if( istep%20==0 )
      std::cout<<istep<<std::endl;
      //print*,minval(h),sum(h)/SIZE,maxval(h)

  }

  std::cout<<"t Step1: Initialize         = "<<std::setw(15)<<Tmr_Step1_Initialize.elapsed()         <<" microseconds"<<std::endl;                 
  std::cout<<"t Step2: DetermineReceivers = "<<std::setw(15)<<Tmr_Step2_DetermineReceivers.elapsed() <<" microseconds"<<std::endl;                         
  std::cout<<"t Step3: DetermineDonors    = "<<std::setw(15)<<Tmr_Step3_DetermineDonors.elapsed()    <<" microseconds"<<std::endl;                      
  std::cout<<"t Step4: GenerateStack      = "<<std::setw(15)<<Tmr_Step4_GenerateStack.elapsed()      <<" microseconds"<<std::endl;                    
  std::cout<<"t Step5: FlowAcc            = "<<std::setw(15)<<Tmr_Step5_FlowAcc.elapsed()            <<" microseconds"<<std::endl;              
  std::cout<<"t Step6: Uplift             = "<<std::setw(15)<<Tmr_Step6_Uplift.elapsed()             <<" microseconds"<<std::endl;             
  std::cout<<"t Step7: Erosion            = "<<std::setw(15)<<Tmr_Step7_Erosion.elapsed()            <<" microseconds"<<std::endl;              
  std::cout<<"t Overall                   = "<<std::setw(15)<<Tmr_Overall.elapsed()                  <<" microseconds"<<std::endl;        


  PrintDEM("out_BW.dem", h, WIDTH, HEIGHT);

  return 0;
}
