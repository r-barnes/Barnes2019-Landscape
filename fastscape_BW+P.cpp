//This file contains a parallel implementation of Braun and Willett's FastScape
//algorithm. The serial implementation was developed by adapting Fortran code
//provided by Braun and attempts to be a faithful reproduction of the ideas
//therein. The parallel version is developed based on ideas expressed in Braun
//and Willett (2013, doi: doi:10.1016/j.geomorph.2012.10.008).
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

  //Each edge node of the DEM is the root of a stack. The `stack_start` array
  //indicates where these nodes are located in the global stack. This is then
  //used to induce parallelism. Note that, especially at the beginning of the
  //program, there may be interior nodes with in depressions or with no local
  //gradient. These may become the starts of stacks as well, so the
  //`stack_start` array may be larger than just the perimeter of the DEM.
  std::vector<int>    stack_start;
  stack_start.reserve(2*WIDTH+2*HEIGHT);
  

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
  for(int y=0;y<HEIGHT;y++)
  for(int x=0;x<WIDTH;x++){
    int c = y*WIDTH+x;
    h[c]  = rand()/(double)RAND_MAX;
    if(x == 0 || y==0 || x==WIDTH-1 || y==HEIGHT-1)
      h[c] = 0;
  }

  //! begining of time stepping
  for(int istep=0;istep<nstep;istep++){
    //The `stack_start` array has an unpredictable size, so we fill it
    //dynamically and reset it each time. This should have minimal impact on the
    //algorithm's speed since the memory is not actually reallocated.
    stack_start.clear();

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
        //Make a note that this is the start of a stack. Later we will use these
        //notes to induce parallelism.
        stack_start.push_back(nstack);
        stack[nstack++] = c;
        find_stack(c,donor,ndon,SIZE,stack,nstack);
      }
    }
    //We add an additional note to the end of `stack_start` that serves as an
    //upper bound on the locations of the cells in the final stack. See below.
    stack_start.push_back(SIZE);

    //! computing drainage area
    for(int i=0;i<SIZE;i++)
      accum[i] = dx*dy;

    for(int s=SIZE-1;s>=0;s--){
      const int c = stack[s];
      if(rec[c]!=c){
        const int n = rec[c];
        accum[n] = accum[n] + accum[c];
      }
    }

    //! adding uplift to landscape
    for(int y=1;y<HEIGHT-1;y++)
    for(int x=1;x<WIDTH-1;x++){
      int c = y*WIDTH+x;
      h[c] += ueq*dt;
    }

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
    #pragma omp parallel        //Spin up the thread team
    #pragma omp single nowait   //The following for loop is executed by only one thread
    for(int ss=0;ss<stack_start.size()-1;ss++){
      //Execute the inner for loop as a task
      const int sstart = stack_start.at(ss);
      const int send   = stack_start.at(ss+1);
      #pragma omp task default(none) shared(stack,rec,accum,length,h) firstprivate(sstart,send)
      for(int s=sstart;s<send;s++){
        const int c = stack[s]; //Cell from which flow originates
        const int n = rec[c];   //Cell receiving the flow
        if(n==c)
          continue;
        const double fact = keq*dt*std::pow(accum[c],meq)/std::pow(length[c],neq);
        const double h0   = h[c];
        double hp         = h0;
        double diff       = 2*tol;
        while(std::abs(diff)>tol){
          h[c] -= (h[c]-h0+fact*std::pow(h[c]-h[n],neq))/(1.+fact*neq*std::pow(h[c]-h[n],neq-1));
          diff  = h[c] - hp;
          hp    = h[c];
        }
      }
    }

    if( istep%20==0 )
      std::cout<<istep<<std::endl;
      //print*,minval(h),sum(h)/SIZE,maxval(h)

  }

  PrintDEM("out_BW+P.dem", h, WIDTH, HEIGHT);

  return 0;
}
