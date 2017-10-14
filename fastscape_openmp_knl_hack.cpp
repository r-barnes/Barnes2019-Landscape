//real, dimension(:), allocatable :: h,accum,length
//real, dimension(:,:), allocatable :: x,y,z
//integer, dimension(:), allocatable :: rec,ndon,stack
//integer, dimension(:,:), allocatable :: donor

//integer WIDTH,HEIGHT,SIZE,nstep,nfreq,nstack
//integer i,j,ij,ii,jj,iii,jjj,ijk,ijr,istep
//real xl,yl,dx,dy,dt,k,n,m,u,l,slope,smax
//real diff,fact,h0,hp,tol

//-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Wno-unused


#include <cstdlib>
#include <limits>
#include <iostream>
#include <cmath>
#include <fstream>
#include <queue>
#include <functional>

const double keq = 2e-6;
const double neq = 2;
const double meq = 0.8;
const double u   = 2e-3;
const double dt  = 1000.;

double *h;
double *accum;
double *length;
int    *rec;
int    *ndon;
int    *stack;
int    *donor;



void PrintDEM(const std::string filename, const int width, const int height){
  std::ofstream fout(filename.c_str());
  fout<<"cols"<<width<<"\n";
  fout<<"nrows"<<height<<"\n";
  fout<<"xllcorner 637500.000\n";
  fout<<"yllcorner 206000.000\n";
  fout<<"cellsize 500.000\n";
  fout<<"NODATA_value -9999\n";
  for(int y=0;y<height;y++){
    for(int x=0;x<width;x++)
      fout<<h[y*width+x]<<" ";
    fout<<"\n";
  }
}



void ErodePoint( int c){
  // #pragma omp critical
  // {
  //   std::cout<<c<<std::endl;
  // }

  //int c=1337;

  if(c!=rec[c]){

    const double fact = keq*dt*std::pow(accum[c],meq)/std::pow(length[c],neq);
    //Use Newton's method to solve backward Euler equation. Fix number of loops
    //to 5, which should be sufficient
    const double hn = h[rec[c]];
    const double h0 = h[c];
    double hnew     = h0;

    for(int i=0;i<5;i++)
      hnew -= (hnew-h0+fact*std::pow(hnew-hn,neq))/(1.+fact*neq*std::pow(hnew-hn,neq-1));

    h[c] = hnew;
  }

  if(ndon[c]>0){
    for(int k=1;k<ndon[c];k++){
      #pragma omp task shared(h)
      ErodePoint(donor[8*c+k]);
    }
    ErodePoint(donor[8*c+0]);
  }
}





int main(){

  //! defining size of the problem
  const int WIDTH    = 501;
  const int HEIGHT   = 501;
  const int SIZE     = WIDTH*HEIGHT;

  //!    allocating memory
  h      = new double[SIZE];
  accum  = new double[SIZE];
  length = new double[SIZE];
  rec    = new int[SIZE];
  ndon   = new int[SIZE];
  stack  = new int[SIZE];
  donor  = new int[8*SIZE];
  
  //double x[SIZE];
  //double y[SIZE];
  //double z[SIZE];

  //Neighbours
  const int nshift[8] = {-1,-WIDTH-1,-WIDTH,-WIDTH+1,1,WIDTH+1,WIDTH,WIDTH-1};
  const double SQRT2 = 1.414213562373095048801688724209698078569671875376948;
  const double dr[8] = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};

  //defining geometrical and temporal constants
  const double xl = 100.e3;
  const double yl = 100.e3;
  const double dx = xl/(WIDTH-1);
  const double dy = yl/(HEIGHT-1);
  //const int nstep = 300;
  const int nstep = 120;

  //! generating initial topography
  for(int y=0;y<HEIGHT;y++)
  for(int x=0;x<WIDTH;x++){
    int ij = y*WIDTH+x;
    h[ij]  = rand()/(double)RAND_MAX;
    if(x == 0 || y==0 || x==WIDTH-1 || y==HEIGHT-1)
      h[ij] = 0;
  }

  //! begining of time stepping
  for(int istep=0;istep<nstep;istep++){

    //! initializing rec and length
    for(int i=0;i<SIZE;i++){
      rec[i]    = i;
      length[i] = 0;
    }

    //! computing receiver array
    for(int y=1;y<HEIGHT-1;y++)
    for(int x=1;x<WIDTH-1;x++){
      int c            = y*HEIGHT+x;
      double slope_max = -std::numeric_limits<double>::infinity(); //Slope max
      int max_n        = c;
      int slope_n      = 0;
      for(int n=0;n<8;n++){
        double slope = (h[c] - h[c+nshift[n]])/dr[n];
        if(slope>slope_max){
          slope_max = slope;
          max_n = n;
          slope_n=dr[n];
        }
      }
      rec[c]    = max_n;
      length[c] = slope_n;
    }

    //! initialising number of donors per node to 0
    for(int i=0;i<SIZE;i++)
      ndon[i] = 0;

    //! computing donor arrays
    for(int i=0;i<SIZE;i++){
      if(rec[i]==i)
        continue;
      const int n        = rec[i];
      donor[8*n+ndon[n]] = i;
      ndon[n]++;
    }

    //! computing drainage area
    for(int i=0;i<SIZE;i++)
      accum[i] = dx*dy;

    for(int s=SIZE-1;s>=0;s--){
      const int c = stack[s];
      if(rec[c]!=c)
        accum[rec[c]] = accum[rec[c]] + accum[c];
    }

    //! adding uplift to landscape
    for(int y=1;y<HEIGHT-1;y++)
    for(int x=1;x<WIDTH-1;x++){
      int c = y*WIDTH+x;
      h[c] += u*dt;
    }

    //! computing stack

    std::queue<int> qu;

    //Initialize stack with cells which do not have dependencies
    #pragma omp parallel
    {
      #pragma omp single nowait
      {
        for(int c=0;c<SIZE;c++){
          if(rec[c]==c){
            #pragma omp task shared(h)
            ErodePoint(c);
          }
        }
      }
    }


    if( istep%20==0 )
      std::cout<<istep<<std::endl;

  }

  PrintDEM("out.dem", WIDTH, HEIGHT);

  return 0;
}
