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

  std::vector<double>       accum;
  std::vector<int>          rec;
  std::vector<unsigned int> ndon;
  std::vector<unsigned int> stack;
  std::vector<unsigned int> donor;
  std::vector<unsigned int> levels;
  
  int nshift[8];

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

  void GenerateStack(){
    unsigned int nstack = 0;
    unsigned int qpoint = 0;

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
      for(unsigned int k=0;k<ndon[c];k++){
        const auto n    = donor[8*c+k];
        stack[nstack++] = n;
      }

      qpoint++;
      if(qpoint==levels.back() && nstack!=levels.back())
        levels.push_back(nstack); //Starting a new level      
    }
    // std::cerr<<"nstack final = "<<nstack<<std::endl;
  }

  void ComputeDraingeArea(){
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

  void AddUplift(){
    //! adding uplift to landscape
    for(int y=1;y<height-1;y++)
    for(int x=1;x<width-1;x++){
      const int c = y*width+x;
      h[c] += ueq*dt;
    }    
  }

  void Erode(){
    auto *hvec     = h.data();
    auto *accvec   = accum.data();
    auto *recvec   = rec.data();
    auto *stackvec = stack.data();

    #pragma acc data default(none) copy(hvec[0:h.size()]) copyin(accvec[0:accum.size()], recvec[0:rec.size()], stackvec[0:stack.size()]) present(dr, nshift) copyin(NO_FLOW)
    for(unsigned int li=0;li<levels.size()-1;li++){
      const auto lvlstart = levels[li];
      const auto lvlend   = levels[li+1];
      #pragma acc parallel loop default(none) //present(hvec,accvec,recvec,dr,nshift)
      for(auto si=lvlstart;si<lvlend;si++){
        const int c = stackvec[si]; //Cell from which flow originates
        if(recvec[c]!=NO_FLOW){  //Can't use continue inside of OpenAcc
          const int n         = c+nshift[recvec[c]];   //Cell receiving the flow
          const double length = dr[recvec[c]];
          const double fact   = keq*dt*std::pow(accvec[c],meq)/std::pow(length,neq);
          const double h0     = hvec[c];
          const double hn     = hvec[n];
//          double hp           = h0;
//          double diff         = 2*tol;
          double hnew         = h0;
          for(int coni=0;coni<10;coni++){
            //Use Newton's method to solve backward Euler equation. Fix number of loops
            //to 5, which should be sufficient
            //for(int i=0;i<5;i++)
            hnew -= (hnew-h0+fact*std::pow(hnew-hn,neq))/(1.+fact*neq*std::pow(hnew-hn,neq-1));
//            diff  = hnew - hp;
//            hp    = hnew;
          }
          hvec[c] = hnew;
        }
      }
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

    nshift[0] = -1;
    nshift[1] = -width-1;
    nshift[2] = -width;
    nshift[3] = -width+1;
    nshift[4] = 1;
    nshift[5] = width+1;
    nshift[6] = width;
    nshift[7] = width-1;

    #pragma acc enter data copyin(dr[0:8], nshift[0:8])
  }

  ~TerrainMorpher(){
    #pragma acc exit data delete(dr, nshift)
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
    for(int istep=0;istep<nstep;istep++){

      ComputeReceivers();

      ComputeDonors();

      GenerateStack();

      ComputeDraingeArea();

      AddUplift();

      // std::cerr<<"Levels: ";
      // for(auto &l: levels)
      //   std::cerr<<l<<" ";
      // std::cerr<<std::endl;

      Erode();

      if( istep%20==0 )
        std::cout<<istep<<std::endl;
        //print*,minval(h),sum(h)/SIZE,maxval(h)
    }
  }
};







int main(){
  //feenableexcept(FE_ALL_EXCEPT);

  const int width  = 501;
  const int height = 501;
  const int nstep  = 120;
  
  TerrainMorpher tm(width,height);
  tm.generateRandomTerrain();


  Timer tmr;
  tm.run(nstep);
  std::cout<<"Calculation time = "<<tmr.elapsed()<<std::endl;

  PrintDEM("out.dem", tm.h, width, height);

  return 0;
}
