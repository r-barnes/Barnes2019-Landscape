//real, dimension(:), allocatable :: h,a,length
//real, dimension(:,:), allocatable :: x,y,z
//integer, dimension(:), allocatable :: rec,ndon,stack
//integer, dimension(:,:), allocatable :: donor

//integer WIDTH,HEIGHT,SIZE,nstep,nfreq,nstack
//integer i,j,ij,ii,jj,iii,jjj,ijk,ijr,istep
//real xl,yl,dx,dy,dt,k,n,m,ueq,l,slope,smax
//real diff,fact,h0,hp,tol


#include <cstdlib>
#include <limits>
#include <iostream>
#include <cmath>
#include <fstream>
#include <fenv.h> //Used to catch floating point NaN issues

constexpr double DINFTY = std::numeric_limits<double>::infinity();

void find_stack(
  const int c, int donor[], int ndon[], int SIZE, int stack[], int &nstack
){
  for(int k=0;k<ndon[c];k++){
    int n           = donor[8*c+k];
    stack[nstack++] = n;
    find_stack(n,donor,ndon,SIZE,stack,nstack);
  }
}

void PrintDEM(const std::string filename, const double h[], const int width, const int height){
  std::ofstream fout(filename.c_str());
  fout<<"ncols "<<width<<"\n";
  fout<<"nrows "<<height<<"\n";
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


int main(){
  unsigned long long cells_processed = 0;
  unsigned long long cells_eroded    = 0;

  //feenableexcept(FE_ALL_EXCEPT);

  //! defining size of the problem
  const int WIDTH = 501;
  const int HEIGHT = 501;
  //const int WIDTH = 101;
  //const int HEIGHT = 101;
  const int SIZE = WIDTH*HEIGHT;

  //!    allocating memory
  double *h      = new double[SIZE];
  double *a      = new double[SIZE];
  double *length = new double[SIZE];
  int    *rec    = new int[SIZE];
  int    *ndon   = new int[SIZE];
  int    *stack  = new int[SIZE];
  int    *donor  = new int[8*SIZE];
  

  //Neighbours
  const int nshift[8] = {-1,-WIDTH-1,-WIDTH,-WIDTH+1,1,WIDTH+1,WIDTH,WIDTH-1};
  const double SQRT2  = 1.414213562373095048801688724209698078569671875376948;
  const double dr[8]  = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};

  //defining geometrical and temporal constants
  const double xl   = 100.e3;
  const double yl   = 100.e3;
  const double dx   = xl/(WIDTH-1);
  const double dy   = yl/(HEIGHT-1);
  //const double dt = 10000.;
  const double dt   = 1000.;
  //const int nstep = 300;
  const int nstep   = 120;
  const double tol  = 1.e-3;

  //! generating initial topography
  for(int y=0;y<HEIGHT;y++)
  for(int x=0;x<WIDTH;x++){
    int c = y*WIDTH+x;
    h[c]  = rand()/(double)RAND_MAX;
    //h[c]=500;
    if(x == 0 || y==0 || x==WIDTH-1 || y==HEIGHT-1)
      h[c] = 0;
  }

  const double k   = 5e-4;
  const double neq = 1;
  const double meq = 0.5;
  const double ueq = 2e-4;

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

    //! initialising number of donors per node to 0
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

    //! computing stack

    int nstack=0;
    for(int c=0;c<SIZE;c++){
      if(rec[c]==c){
        stack[nstack++] = c;
        find_stack(c,donor,ndon,SIZE,stack,nstack);
      }
    }

    //! computing drainage area
    for(int i=0;i<SIZE;i++)
      a[i] = dx*dy;

    for(int s=SIZE-1;s>=0;s--){
      const int c = stack[s];
      if(rec[c]!=c)
        a[rec[c]] = a[rec[c]] + a[c];
    }

    //! adding uplift to landscape
    for(int y=1;y<HEIGHT-1;y++)
    for(int x=1;x<WIDTH-1;x++){
      int c = y*WIDTH+x;
      h[c] += ueq*dt;
    }

    std::cout<<"rec: ";
    for(int c=3*WIDTH;c<4*WIDTH;c++){
      if(rec[c]==c)
        std::cout<<-1<<" ";
      else
        std::cout<<rec[c]<<" ";
    }
    std::cout<<std::endl;

    //! computing erosion
    std::cerr<<"length: ";
    for(int i=3*WIDTH;i<4*WIDTH;i++)
      std::cerr<<length[i]<<" ";
    std::cerr<<std::endl;

    for(int s=0;s<SIZE;s++){
      cells_processed++;
      const int c = stack[s]; //Cell from which flow originates
      const int n = rec[c];   //Cell receiving the flow
      if(n==c)
        continue;
      cells_eroded++;
      const double fact = k*dt*std::pow(a[c],meq)/std::pow(length[c],neq);
      const double h0   = h[c];
      double hp         = h0;
      double diff       = 2*tol;
      while(std::abs(diff)>tol){
        //Use Newton's method to solve backward Euler equation. Fix number of loops
        //to 5, which should be sufficient
        //for(int i=0;i<5;i++)
        h[c] -= (h[c]-h0+fact*std::pow(h[c]-h[n],neq))/(1.+fact*neq*std::pow(h[c]-h[n],neq-1));
        diff  = h[c] - hp;
        hp    = h[c];
      }
    }

    if( istep%20==0 )
      std::cout<<istep<<std::endl;
      //print*,minval(h),sum(h)/SIZE,maxval(h)

  }

  std::cout<<"Cells processed = "<<cells_processed<<std::endl;
  std::cout<<"Cells eroded = "<<cells_eroded<<std::endl;

  PrintDEM("out.dem", h, WIDTH, HEIGHT);

  return 0;
}
