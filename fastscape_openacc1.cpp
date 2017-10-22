#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <array>

const double keq     = 2e-6;
const double neq     = 2;
const double meq     = 0.8;
const double ueq     = 2e-3;
const double dt      = 1000.;
const int    NO_FLOW = -1;

constexpr double DINFTY  = std::numeric_limits<double>::infinity();

const double SQRT2  = 1.414213562373095048801688724209698078569671875376948;
const double dr[8]  = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};



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





void GenerateRandomTerrain(double *const h, const int width, const int height){
  for(int y=0;y<height;y++)
  for(int x=0;x<width;x++){
    const int c = y*width+x;
    h[c]  = rand()/(double)RAND_MAX;
    if(x == 0 || y==0 || x==width-1 || y==height-1)
      h[c] = 0;
  }
}



void ComputeReceivers(
  const double *const h,
  const int    *const nshift,
  const int width,
  const int height,
        int    *const rec
){
  const int size = width*height;

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


void ComputeDonors(
  const int *const rec,
  const int *const nshift,
  const int width,
  const int height,
  int *const donor,
  int *const ndon
){
  const int size = width*height;

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



void GenerateStack(
  const int *const rec,
  const int *const donor,
  const int *const ndon,
  const int width,
  const int height,
  std::vector<int> &levels,
  int *const stack
){
  const int size = width*height;

  int nstack = 0;
  int qpoint = 0;

  levels.clear();
  levels.push_back(0);

  //Load cells without dependencies into the queue
  for(int c=0;c<size;c++){
    if(rec[c]==NO_FLOW)
      stack[nstack++] = c;
  }
  levels.push_back(nstack); //Last cell of this level

  while(qpoint<nstack){
    const auto c = stack[qpoint];
    for(int k=0;k<ndon[c];k++){
      const auto n    = donor[8*c+k];
      stack[nstack++] = n;
    }

    qpoint++;
    if(qpoint==levels.back() && nstack!=levels.back())
      levels.push_back(nstack); //Starting a new level      
  }
  // std::cerr<<"nstack final = "<<nstack<<std::endl;
}



void ComputeDraingeArea(
  const int *const rec,
  const int *const nshift,
  const int *const stack,
  const int width,
  const int height,
  double *const accum
){
  const int size = width*height;

  const double xl   = 100.e3;
  const double yl   = 100.e3;
  const double dx   = xl/(width-1);
  const double dy   = yl/(height-1);
  const double area = dx*dy;

  //! computing drainage area
  for(int i=0;i<size;i++)
    accum[i] = area;

  for(int s=size-1;s>=0;s--){
    const int c = stack[s];
    if(rec[c]!=NO_FLOW){
      const int n = c+nshift[rec[c]];
      accum[n]   += accum[c];
    }
  }    
}



void AddUplift(
  const int width,
  const int height,
  double *const h
){
  //! adding uplift to landscape
  for(int y=1;y<height-1;y++)
  for(int x=1;x<width-1;x++){
    const int c = y*width+x;
    h[c] += ueq*dt;
  }    
}



void Erode(
  const double *const accum,
  const int    *const rec,
  const int    *const stack,
  const int    *const nshift,
  const std::vector<int> levels,
  const int width,
  const int height,
  double *const h
){
  const int size = width*height;

  const double tol = 1.e-3;

  #pragma acc data default(none) copy(hvec[0:h.size()]) copyin(accvec[0:accum.size()], recvec[0:rec.size()], stackvec[0:stack.size()]) present(dr, nshift)
  for(unsigned int li=0;li<levels.size()-1;li++){
    const auto lvlstart = levels[li];
    const auto lvlend   = levels[li+1];
    #pragma acc parallel loop default(none) //present(hvec,accvec,recvec,dr,nshift)
    for(auto si=lvlstart;si<lvlend;si++){
      const int c = stack[si]; //Cell from which flow originates
      if(rec[c]!=NO_FLOW){  //Can't use continue inside of OpenAcc
        const int n         = c+nshift[rec[c]];   //Cell receiving the flow
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
  }
}



double* TerrainMorpher(
  const int width,
  const int height,
  const int nstep
){
  const int size = width*height;

  double *const h     = new double[size];
  double *const accum = new double[size];
  int    *const rec   = new int[size];
  int    *const ndon  = new int[size];
  int    *const stack = new int[size];
  int    *const donor = new int[8*size];
  std::vector<int> levels;

  const int nshift[8] = {-1,-width-1,-width,-width+1,1,width+1,width,width-1};
  
  GenerateRandomTerrain(h, width, height);

  #pragma acc enter data copyin(dr[0:8], nshift[0:8])

  for(int step=0;step<nstep;step++){

    ComputeReceivers(h, nshift, width, height, rec);

    ComputeDonors(rec, nshift, width, height, donor, ndon);

    GenerateStack(rec, donor, ndon, width, height, levels, stack);

    ComputeDraingeArea(rec, nshift, stack, width, height, accum);

    AddUplift(width, height, h);

    Erode(accum, rec, stack, nshift, levels, width, height, h);

    if( step%20==0 )
      std::cout<<step<<std::endl;
  }


  #pragma acc exit data delete(dr,nshift)

  delete[] accum;
  delete[] rec;
  delete[] ndon;
  delete[] stack;
  delete[] donor;

  return h;
}



int main(){
  //feenableexcept(FE_ALL_EXCEPT);

  const int width  = 501;
  const int height = 501;
  const int nstep  = 120;

  auto h = TerrainMorpher(width,height,nstep);

  PrintDEM("out.dem", h, width, height);

  delete[] h;

  return 0;
}
