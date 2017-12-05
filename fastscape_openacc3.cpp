#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include <array>
#include <random>
#include "Timer.hpp"

const double keq     = 2e-6;
const double neq     = 2;
const double meq     = 0.8;
const double ueq     = 2e-3;
const double dt      = 1000.;
const int    NO_FLOW = -1;

constexpr double DINFTY  = std::numeric_limits<double>::infinity();

const double SQRT2  = 1.414213562373095048801688724209698078569671875376948;
const double dr[8]  = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2};


#include "draw_stuff.hpp"


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
  srand(std::random_device()());
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
  const int           width,
  const int           height,
        int    *const rec
){
  const int size = width*height;

  #pragma acc kernels default(none) present(rec,nshift,dr,h)
  {
    //! initializing rec and length
    for(int i=0;i<size;i++)
      rec[i] = NO_FLOW;

    //! computing receiver array
    #pragma acc loop collapse(2) independent
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
      const int c      = y*width+x;
      double max_slope = -DINFTY; //TODO: Wrong, in the case of flats
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

  #pragma acc kernels default(none) present(ndon,donor,rec,nshift)
  {
    //! initialising number of donors per node to 0
    for(int i=0;i<size;i++)
      ndon[i] = 0;

    //! computing donor arrays
    #pragma acc loop independent
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW)
        continue;
      const int n        = c+nshift[rec[c]];
      int this_n;
      #pragma acc atomic capture
      this_n = ndon[n]++;
      donor[8*n+this_n] = c;
    }    
  }
}



void GenerateStack(
  const int *const rec,
  const int *const donor,
  const int *const ndon,
  const int width,
  const int height,
  int *const levels,
  int &lvlsize,
  int *const stack
){
  const int size = width*height;

  int nstack = 0;
  int qpoint = 0;

  #pragma acc kernels default(none) present(levels,stack,ndon,donor,rec)
  {

    levels[0] = 0;
    lvlsize   = 1;

    //Load cells without dependencies into the queue
    for(int c=0;c<size;c++){
      if(rec[c]==NO_FLOW)
        stack[nstack++] = c;
    }
    levels[lvlsize++] = nstack; //Last cell of this level

    while(qpoint<nstack){
      const auto c = stack[qpoint];
      for(int k=0;k<ndon[c];k++){
        const auto n    = donor[8*c+k];
        stack[nstack++] = n;
      }

      //TODO: What's this about, then?
      qpoint++;
      if(qpoint==levels[lvlsize-1] && nstack!=levels[lvlsize-1])
        levels[lvlsize++] = nstack; //Starting a new level      
    }
    // std::cerr<<"nstack final = "<<nstack<<std::endl;
  }
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

  //#pragma acc kernels default(none) present(accum,nshift,rec,stack)

  #pragma acc update host(rec,stack)

  {
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

  #pragma acc update device(accum)
}



void AddUplift(
  const int width,
  const int height,
  double *const h
){
  const int size = width*height;

  //! adding uplift to landscape
  #pragma acc kernels default(none) async present(h)
  #pragma acc loop collapse(2) independent
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
  const int    *const levels,
  const int lvlsize,
  const int width,
  const int height,
  double *const h
){
  const int size = width*height;

  const double tol = 1.e-3;

  #pragma acc wait

  //TODO: I think we can skip all of level 0!

  //#pragma acc data default(none) copyin(accum[0:size], rec[0:size], stack[0:size]) present(dr, nshift, h)
  #pragma acc data default(none) present(dr, nshift, h, accum, rec, stack)
  for(unsigned int li=0;li<lvlsize-1;li++){
    const auto lvlstart = levels[li];
    const auto lvlend   = levels[li+1];
    #pragma acc parallel loop default(none) //present(hvec,accvec,recvec,dr,nshift)
    for(auto si=lvlstart;si<lvlend;si++){
      const int c = stack[si]; //Cell from which flow originates
      if(rec[c]!=NO_FLOW){  //Can't use continue inside of OpenAcc
        const int n         = c+nshift[rec[c]];   //Cell receiving the flow

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
}



double* TerrainMorpher(
  const int width,
  const int height,
  const int nstep
){
  const int size = width*height;

  double *const h      = new double[size];
  double *const accum  = new double[size];
  int    *const rec    = new int[size];
  int    *const ndon   = new int[size];
  int    *const stack  = new int[size];
  int    *const donor  = new int[8*size];
  int    *const levels = new int[2*width]; //Pray this is enough
  int    lvlsize;

  const int nshift[8] = {-1,-width-1,-width,-width+1,1,width+1,width,width-1};
  
  GenerateRandomTerrain(h, width, height);

  #pragma acc enter data create(rec[0:size], stack[0:size], accum[0:size], ndon[0:size], donor[0:8*size], levels[0:2*width]) copyin(dr[0:8], nshift[0:8], h[0:size])

  for(int step=0;step<=nstep;step++){

    ComputeReceivers(h, nshift, width, height, rec);

    AddUplift(width, height, h);

    ComputeDonors(rec, nshift, width, height, donor, ndon);

    GenerateStack(rec, donor, ndon, width, height, levels, lvlsize, stack);

    ComputeDraingeArea(rec, nshift, stack, width, height, accum);

    Erode(accum, rec, stack, nshift, levels, lvlsize, width, height, h);

    if( step%20==0 )
      std::cout<<step<<std::endl;
  }

  #pragma acc exit data copyout(h[0:size]) delete(dr,nshift,rec,stack,accum)

  StartPicture();
  DrawQueue(h, width, height, nshift, rec, stack, levels, lvlsize);
  EndPicture();

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

  Timer tmr;
  auto h = TerrainMorpher(width,height,nstep);
  std::cout<<"Calculation time = "<<tmr.elapsed()<<std::endl;

  //PrintDEM("out_openacc3.dem", h, width, height);

  delete[] h;

  return 0;
}
