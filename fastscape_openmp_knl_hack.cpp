#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

const double keq = 2e-6;
const double neq = 2;
const double meq = 0.8;
const double ueq = 2e-3;
const double dt  = 1000.;

constexpr double DINFTY  = std::numeric_limits<double>::infinity();
const     int    NO_FLOW = -1;

const double SQRT2  = 1.414213562373095048801688724209698078569671875376948;
const double dr[8]  = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};



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



void find_stack(
  const int c, 
  const std::vector<int> &donor,
  const std::vector<int> &ndon,
  std::vector<int> &stack,
  int &nstack
){
  for(int k=0;k<ndon[c];k++){
    int n           = donor[8*c+k];
    stack[nstack++] = n;
    find_stack(n,donor,ndon,stack,nstack);
  }
}



void ErodePoint(
  const int c,
  const std::vector<int>    &donor,
  const std::vector<int>    &ndon,
  const std::vector<int>    &rec,
  const std::vector<double> &accum,
  const std::vector<int>    &nshift,
        std::vector<double> &h
){
  // #pragma omp critical
  // {
  //   std::cout<<c<<std::endl;
  // }

  if(rec[c]!=NO_FLOW){
    const int    n      = c+nshift[rec[n]];
    const double length = dr[rec[n]];
    const double fact   = keq*dt*std::pow(accum[c],meq)/std::pow(length,neq);
    //Use Newton's method to solve backward Euler equation. Fix number of loops
    //to 5, which should be sufficient
    const double hn = h[n];
    const double h0 = h[c];
    double hnew     = h0;
    for(int i=0;i<5;i++)
      hnew -= (hnew-h0+fact*std::pow(hnew-hn,neq))/(1.+fact*neq*std::pow(hnew-hn,neq-1));

    h[c] = hnew;
  }

  if(ndon[c]>0){
    for(int k=1;k<ndon[c];k++){
      #pragma omp task shared(h)
      ErodePoint(donor[8*c+k],donor,ndon,rec,accum,nshift,h);
    }
    ErodePoint(donor[8*c+0],donor,ndon,rec,accum,nshift,h);
  }
}





int main(){

  //! defining size of the problem
  const int WIDTH    = 501;
  const int HEIGHT   = 501;
  const int SIZE     = WIDTH*HEIGHT;

  //!    allocating memory
  std::vector<double> h    (  SIZE);
  std::vector<double> accum(  SIZE);
  std::vector<int>    rec  (  SIZE);
  std::vector<int>    ndon (  SIZE);
  std::vector<int>    stack(  SIZE);
  std::vector<int>    donor(8*SIZE);
  

  //Neighbours
  const std::vector<int> nshift= {{-1,-WIDTH-1,-WIDTH,-WIDTH+1,1,WIDTH+1,WIDTH,WIDTH-1}};

  //defining geometrical and temporal constants
  const double xl = 100.e3;
  const double yl = 100.e3;
  const double dx = xl/(WIDTH-1);
  const double dy = yl/(HEIGHT-1);
  const int nstep = 120;

  //! generating initial topography
  for(int y=0;y<HEIGHT;y++)
  for(int x=0;x<WIDTH;x++){
    const int c = y*WIDTH+x;
    h[c]  = rand()/(double)RAND_MAX;
    if(x == 0 || y==0 || x==WIDTH-1 || y==HEIGHT-1)
      h[c] = 0;
  }

  //! begining of time stepping
  for(int istep=0;istep<nstep;istep++){

    //! initializing rec and length
    for(int i=0;i<SIZE;i++)
      rec[i] = NO_FLOW;

    //! computing receiver array
    for(int y=1;y<HEIGHT-1;y++)
    for(int x=1;x<WIDTH-1;x++){
      const int c      = y*WIDTH+x;
      double max_slope = -DINFTY;
      int    max_n     = NO_FLOW;
      for(int n=0;n<8;n++){
        double slope = (h[c] - h[c+nshift[n]])/dr[n];
        if(slope>max_slope){
          max_slope = slope;
          max_n = n;
        }
      }
      rec[c]    = max_n;
    }

    //! initialising number of donors per node to 0
    for(int i=0;i<SIZE;i++)
      ndon[i] = 0;

    //! computing donor arrays
    for(int c=0;c<SIZE;c++){
      if(rec[c]==NO_FLOW)
        continue;
      const int n        = c+nshift[rec[c]];
      donor[8*n+ndon[n]] = c;
      ndon[n]++;
    }

    //! computing stack

    int nstack=0;
    for(int c=0;c<SIZE;c++){
      if(rec[c]==NO_FLOW){
        stack[nstack++] = c;
        find_stack(c,donor,ndon,stack,nstack);
      }
    }

    //! computing drainage area
    for(int i=0;i<SIZE;i++)
      accum[i] = dx*dy;

    for(int s=SIZE-1;s>=0;s--){
      const int c = stack[s];
      if(rec[c]!=NO_FLOW){
        const int n = c+nshift[rec[c]];
        accum[n]   += accum[c];
      }
    }

    //! adding uplift to landscape
    for(int y=1;y<HEIGHT-1;y++)
    for(int x=1;x<WIDTH-1;x++){
      int c = y*WIDTH+x;
      h[c] += ueq*dt;
    }



    //Initialize stack with cells which do not have dependencies
    #pragma omp parallel
    {
      #pragma omp single nowait
      {
        for(int c=0;c<SIZE;c++){
          if(rec[c]==NO_FLOW){
            #pragma omp task shared(h)
            ErodePoint(c,donor,ndon,rec,accum,nshift,h);
          }
        }
      }
    }



    if( istep%20==0 )
      std::cout<<istep<<std::endl;
      //print*,minval(h),sum(h)/SIZE,maxval(h)

  }

  PrintDEM("out.dem", h, WIDTH, HEIGHT);

  return 0;
}
