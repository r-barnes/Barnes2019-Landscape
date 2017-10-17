#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <array>

const double keq = 2e-6;
const double neq = 2;
const double meq = 0.8;
const double ueq = 2e-3;
const double dt  = 1000.;

constexpr double DINFTY  = std::numeric_limits<double>::infinity();

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



class TerrainMorpher {
 public:
  std::vector<double> h;

 private:
  static const int NO_FLOW = -1;

  int width;
  int height;
  int size;

  std::vector<double> accum;
  std::vector<int>    rec;
  std::vector<int>    ndon;
  std::vector<int>    stack;
  std::vector<int>    donor;
  
  std::array<int,8> nshift;

  void ComputeReceivers(){
    //! initializing rec and length
    for(int i=0;i<size;i++)
      rec[i] = NO_FLOW;

    //! computing receiver array
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
      const int c      = y*width+x;
      double max_slope = -DINFTY;
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
      const int n     = donor[8*c+k];
      stack[nstack++] = n;
      FindStack(n,nstack);
    }
  }

  void GenerateStack(){
    //! computing stack
    int nstack=0;
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW){
        stack[nstack++] = c;
        FindStack(c, nstack);
      }
    }    
  }

  void ComputeDraingeArea(){
    const double xl   = 100.e3;
    const double yl   = 100.e3;
    const double dx   = xl/(width-1);
    const double dy   = yl/(height-1);
    const double area = dx*dy;

    //! computing drainage area
    for(int i=0;i<size;i++)
      accum[i] = dx*dy;

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

  void Erode(){
    for(int s=0;s<size;s++){
      const int c = stack[s]; //Cell from which flow originates
      if(rec[c]==NO_FLOW)
        continue;
      const int n = c+nshift[rec[c]];   //Cell receiving the flow
      const double length = dr[rec[c]];
      const double fact   = keq*dt*std::pow(accum[c],meq)/std::pow(length,neq);
      const double h0     = h[c];
      const double hn     = h[n];
      double hp           = h0;
      double diff         = 2*tol;
      double hnew         = h0;
      while(std::abs(diff)>tol){
        //Use Newton's method to solve backward Euler equation. Fix number of loops
        //to 5, which should be sufficient
        //for(int i=0;i<5;i++)
        hnew -= (hnew-h0+fact*std::pow(hnew-hn,neq))/(1.+fact*neq*std::pow(hnew-hn,neq-1));
        diff  = hnew - hp;
        hp    = hnew;
      }
      h[c] = hnew;
    }    
  }

 public:
  const double tol = 1.e-3;

 public:
  TerrainMorpher(const int width0, const int height0){
    width  = width0;
    height = height0;
    size   = width*height;

    h.resize    (  size);
    accum.resize(  size);
    rec.resize  (  size);
    ndon.resize (  size);
    stack.resize(  size);
    donor.resize(8*size);

    nshift = {{-1,-width-1,-width,-width+1,1,width+1,width,width-1}};
  }

  void generateRandomTerrain(){
    for(int y=0;y<height;y++)
    for(int x=0;x<width;x++){
      const int c = y*width+x;
      h[c]  = rand()/(double)RAND_MAX;
      if(x == 0 || y==0 || x==width-1 || y==height-1)
        h[c] = 0;
    }
  }

  void run(const int nstep){
    //! begining of time stepping
    for(int istep=0;istep<nstep;istep++){

      ComputeReceivers();

      ComputeDonors();

      GenerateStack();

      ComputeDraingeArea();

      AddUplift();

      Erode();

      if( istep%20==0 )
        std::cout<<istep<<std::endl;
        //print*,minval(h),sum(h)/SIZE,maxval(h)
    }
  }
};







int main(){
  //feenableexcept(FE_ALL_EXCEPT);

  const int width = 501;
  const int height = 501;
  
  auto tm = TerrainMorpher(width,height);
  tm.generateRandomTerrain();

  const int nstep = 120;
  tm.run(nstep);

  PrintDEM("out.dem", tm.h, width, height);

  return 0;
}
