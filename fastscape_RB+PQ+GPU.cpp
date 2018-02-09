#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fenv.h> //Used to catch floating point NaN issues
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include "random.hpp"
#include <vector>
#include "CumulativeTimer.hpp"

#ifndef _OPENMP
  #define omp_get_thread_num()  0
  #define omp_get_num_threads() 1
  #define omp_get_max_threads() 1
#endif

typedef double TerrainType;

void PrintDEM(
  const std::string filename, 
  const double *const h,
  const int width,
  const int height
){
  std::ofstream fout(filename.c_str());
  fout<<"ncols "<<(width- 2)<<"\n";
  fout<<"nrows "<<(height-2)<<"\n";
  fout<<"xllcorner 637500.000\n"; //Arbitrarily chosen value
  fout<<"yllcorner 206000.000\n"; //Arbitrarily chosen value
  fout<<"cellsize 500.000\n";     //Arbitrarily chosen value
  fout<<"NODATA_value -9999\n";
  for(int y=1;y<height-1;y++){
    for(int x=1;x<width-1;x++)
      fout<<h[y*width+x]<<" ";
    fout<<"\n";
  }
}

void GenerateRandomTerrain(TerrainType *const h, const int width, const int height){
  //srand(std::random_device()());
  for(int y=0;y<height;y++)
  for(int x=0;x<width;x++){
    const int c = y*width+x;
    h[c]  = rand()/(double)RAND_MAX;
    if(x == 0 || y==0 || x==width-1 || y==height-1)
      h[c] = 0;
    if(x == 1 || y==1 || x==width-2 || y==height-2)
      h[c] = 0;
  }
}  





TerrainType* FastScape(const int width, const int height, const int nstep){
  const int    NO_FLOW = -1;
  const double SQRT2   = 1.414213562373095048801688724209698078569671875376948;

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

  const int size = width*height;  //Size of DEM (width*height)

  //nshift offsets:
  //1 2 3
  //0   4
  //7 6 5
  const int nshift[8] = {-1,-width-1,-width,-width+1,1,width+1,width,width-1}; //Offset from a focal cell's index to its neighbours

  // CumulativeTimer Tmr_Step1_Initialize;
  // CumulativeTimer Tmr_Step2_DetermineReceivers;
  // CumulativeTimer Tmr_Step3_DetermineDonors;
  // CumulativeTimer Tmr_Step4_GenerateOrder;
  // CumulativeTimer Tmr_Step5_FlowAcc;
  // CumulativeTimer Tmr_Step6_Uplift;
  // CumulativeTimer Tmr_Step7_Erosion;
  // CumulativeTimer Tmr_Overall;


  double *const h     = new double[size]; //Digital elevation model (height)
  double *const accum = new double[size]; //Flow accumulation at each point         
  int    *const rec   = new int[size];    //Index of receiving cell      
  int    *const ndon  = new int[size];    //Indices of a cell's donor cells      
  int    *const donor = new int[8*size];  //How many donors a cell has        

  const int stack_size = 10*size;
  const int level_size = size;
  const int gangs      = 2;

  //TODO: Make smaller, explain max
  const int stack_width = std::max(100,stack_size/gangs);   //Number of stack entries available to each thread
  const int level_width = std::max(100,level_size/gangs);   //Number of level entries available to each thread

  GenerateRandomTerrain(h, width, height);

  int *const stack  = new int[stack_size];

  //It's difficult to know how much memory should be allocated for levels. For
  //a square DEM with isotropic dispersion this is approximately sqrt(E/2). A
  //diagonally tilted surface with isotropic dispersion may have sqrt(E)
  //levels. A tortorously sinuous river may have up to E*E levels. We
  //compromise and choose a number of levels equal to the perimiter because
  //why not?  
  int *const levels = new int[level_size];

  // Tmr_Step1_Initialize.stop();

  #pragma omp target teams num_teams(gangs) default(none) map(alloc:accum[0:size],rec[0:size],ndon[0:size],donor[0:8*size],stack[0:stack_size],levels[0:level_size]) map(tofrom:h[0:size]) map(to:nshift[0:8],dr[0:8]) map(to:width,height,stack_width,level_width) shared(nshift,dr)
  {
    int  nlevel = 0;

    //! initializing rec
    #pragma omp distribute parallel for default(none)
    for(int i=0;i<size;i++)
      rec[i] = NO_FLOW;

    #pragma omp distribute parallel for default(none)
    for(int i=0;i<size;i++)
      ndon[i] = 0;

    for(int step=0;step<=nstep;step++){

      // std::cerr<<"Step = "<<step<<std::endl;

      /////////////////////////////
      //Compute receivers
      #pragma omp distribute parallel for collapse(2) schedule(static) 
      for(int y=2;y<height-2;y++)
      for(int x=2;x<width-2;x++){
        const int c      = y*width+x;

      //   #pragma omp critical
      //   {
      //   std::cout<<"height = "<<height<<std::endl;
      //   std::cout<<"width = "<<width<<std::endl;
      //   std::cout<<"c = "<<c<<std::endl;
      //   std::cout<<"size = "<<size<<std::endl;
      // }

        //The slope must be greater than zero for there to be downhill flow;
        //otherwise, the cell is marekd NO_FLOW
        double max_slope = 0;
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


      /////////////////////////////
      //Compute donors

      //The B&W method of developing the donor array has each focal cell F inform
      //its receiving cell R that F is a donor of R. Unfortunately, parallelizing
      //this is difficult because more than one cell might be informing R at any
      //given time. Atomics are a solution, but they impose a performance cost
      //(though using the latest and greatest hardware decreases this penalty).

      //Instead, we invert the operation. Each focal cell now examines its
      //neighbours to see if it receives from them. Each focal cell is then
      //guaranteed to have sole write-access to its location in the donor array.

      #pragma omp distribute parallel for collapse(2) schedule(static)
      for(int y=1;y<height-1;y++)
      for(int x=1;x<width-1;x++){
        const int c = y*width+x;
        ndon[c] = 0;
        for(int ni=0;ni<8;ni++){
          const int n = c+nshift[ni];
          if(rec[n]!=NO_FLOW && n+nshift[rec[n]]==c){
            donor[8*c+ndon[c]] = n;
            ndon[c]++;
          }
        }
      }

      /////////////////////////////
      //Generate order

      int nstack = 0;

      int *mystack  = &stack [omp_get_team_num()*stack_width];
      int *mylevels = &levels[omp_get_team_num()*level_width];

      // std::cerr<<omp_get_team_num()<<" stack ="<<mystack<<std::endl;
      // std::cerr<<omp_get_team_num()<<" levels ="<<mylevels<<std::endl;


      mylevels[0] = 0;
      nlevel      = 1;

      //Outer edge
      #pragma omp distribute
      for(int y=1;y<height-1;y++){
        mystack[nstack++] = y*width+1;          //assert(nstack<stack_width);
        mystack[nstack++] = y*width+(width-2);  //assert(nstack<stack_width);
      }

      #pragma omp distribute
      for(int x=1;x<width-1;x++){
        mystack[nstack++] =          1*width+x; //assert(nstack<stack_width);
        mystack[nstack++] = (height-2)*width+x; //assert(nstack<stack_width);
      }

      //End of outer edge
      mylevels[nlevel++] = nstack; //Last cell of this level

      //Interior cells
      //TODO: Outside edge is always NO_FLOW. Maybe this can get loaded once?
      //Load cells without dependencies into the queue
      //TODO: Why can't I use nowait here?
      #pragma omp distribute collapse(2)
      for(int y=2;y<height-2;y++)
      for(int x=2;x<width -2;x++){
        const int c = y*width+x;
        if(rec[c]==NO_FLOW){
          mystack[nstack++] = c;      
          // if(nstack>=stack_width){
          //   #pragma omp critical
          //   std::cout<<"nstack = "<<nstack<<" stack_width="<<stack_width<<std::endl;
          // }          
          //assert(nstack<stack_width);
        }
      }
      //Last cell of this level
      mylevels[nlevel++] = nstack;              
      //assert(nlevel<level_width); 

      int level_bottom = -1;
      int level_top    = 0;

      while(level_bottom<level_top){
        level_bottom = level_top;
        level_top    = nstack;
        //#pragma omp parallel for
        for(int si=level_bottom;si<level_top;si++){
          const auto c = mystack[si];
          for(int k=0;k<ndon[c];k++){
            const auto n    = donor[8*c+k];
            mystack[nstack++] = n;         
            // if(nstack>=stack_width)
            //        #pragma omp critical
            //   std::cout<<"nstack = "<<nstack<<" stack_width="<<stack_width<<std::endl;
            //assert(nstack<=stack_width);
          }
        }

        mylevels[nlevel++] = nstack; //Starting a new level      
      }

      //End condition for the loop places two identical entries
      //at the end of the stack. Remove one.
      nlevel--;

    // #pragma omp critical
    // {
    //   for(int i=0;i<nlevel;i++)
    //     std::cerr<<"team "<<omp_get_team_num()<<" levels "<<mylevels[i]<<"-"<<mylevels[i+1]<<std::endl;
    // }

    ////////////////////////////////////////////////
    //Flow Accumulation

    #pragma omp parallel for
    for(int i=mylevels[0];i<mylevels[nlevel-1];i++){
      const int c = mystack[i];
      accum[c] = cell_area;
    }

    for(int li=nlevel-2;li>=0;li--){
      const int lvlstart = mylevels[li];
      const int lvlend   = mylevels[li+1];
      #pragma omp parallel for
      for(int si=lvlstart;si<lvlend;si++){
        const int c = mystack[si];
        if(rec[c]!=NO_FLOW){
          const int n = c+nshift[rec[c]];
          accum[n]   += accum[c];
        }
      }
    }    



    ////////////////////////////////////
    //UPLIFT

    //Start at levels[1] so we don't elevate the outer edge
    #pragma omp parallel for
    for(int i=mylevels[1];i<mylevels[nlevel-1];i++){
      const int c = mystack[i];
      h[c]       += ueq*dt; 
    }

    // std::cerr<<"Eroding..."<<step<<std::endl;


    //#pragma omp parallel default(none)
    for(int li=0;li<nlevel-1;li++){
      const int lvlstart = mylevels[li];
      const int lvlend   = mylevels[li+1];
      // #pragma omp critical
      // std::cerr<<omp_get_team_num()<<"-"<<omp_get_thread_num()<<" going from lvlstart="<<lvlstart<<" to "<<lvlend<<std::endl;
      #pragma omp parallel for default(none) shared(mystack) firstprivate(nshift,dr)
      for(int si=lvlstart;si<lvlend;si++){
        const int c = mystack[si];          //Cell from which flow originates
        if(rec[c]==NO_FLOW)
          continue;
        const int n = c+nshift[rec[c]];    //Cell receiving the flow

        const double length = dr[rec[c]];
        const double fact   = keq*dt*std::pow(accum[c],meq)/std::pow(length,neq);
        const double h0     = h[c];        //Elevation of focal cell
        const double hn     = h[n];        //Elevation of neighbouring (receiving, lower) cell
        double hnew         = h0;          //Current updated value of focal cell
        double hp           = h0;          //Previous updated value of focal cell
        double diff         = 2*tol;       //Difference between current and previous updated values
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


    // Tmr_Step2_DetermineReceivers.start ();   ComputeReceivers  ();                      Tmr_Step2_DetermineReceivers.stop ();
    // Tmr_Step3_DetermineDonors.start    ();   ComputeDonors     ();                      Tmr_Step3_DetermineDonors.stop    ();
    // Tmr_Step4_GenerateOrder.start      ();   GenerateOrder     (stack,levels,nlevel);   Tmr_Step4_GenerateOrder.stop      ();
    // Tmr_Step5_FlowAcc.start            ();   ComputeFlowAcc    (stack,levels,nlevel);   Tmr_Step5_FlowAcc.stop            ();
    // Tmr_Step6_Uplift.start             ();   AddUplift         (stack,levels,nlevel);   Tmr_Step6_Uplift.stop             ();
    // Tmr_Step7_Erosion.start            ();   Erode             (stack,levels,nlevel);   Tmr_Step7_Erosion.stop            ();

    // Tmr_Overall.stop();

    // std::cout<<"t Step1: Initialize         = "<<std::setw(15)<<Tmr_Step1_Initialize.elapsed()         <<" microseconds"<<std::endl;                 
    // std::cout<<"t Step2: DetermineReceivers = "<<std::setw(15)<<Tmr_Step2_DetermineReceivers.elapsed() <<" microseconds"<<std::endl;                         
    // std::cout<<"t Step3: DetermineDonors    = "<<std::setw(15)<<Tmr_Step3_DetermineDonors.elapsed()    <<" microseconds"<<std::endl;                      
    // std::cout<<"t Step4: GenerateOrder      = "<<std::setw(15)<<Tmr_Step4_GenerateOrder.elapsed()      <<" microseconds"<<std::endl;                    
    // std::cout<<"t Step5: FlowAcc            = "<<std::setw(15)<<Tmr_Step5_FlowAcc.elapsed()            <<" microseconds"<<std::endl;              
    // std::cout<<"t Step6: Uplift             = "<<std::setw(15)<<Tmr_Step6_Uplift.elapsed()             <<" microseconds"<<std::endl;             
    // std::cout<<"t Step7: Erosion            = "<<std::setw(15)<<Tmr_Step7_Erosion.elapsed()            <<" microseconds"<<std::endl;              
    // std::cout<<"t Overall                   = "<<std::setw(15)<<Tmr_Overall.elapsed()                  <<" microseconds"<<std::endl;        

    delete[] accum;
    delete[] rec;
    delete[] ndon;
    delete[] donor;
    delete[] stack;
    delete[] levels;

    return h;
  }

int main(int argc, char **argv){
  //feenableexcept(FE_ALL_EXCEPT);

  if(argc!=5){
    std::cerr<<"Syntax: "<<argv[0]<<" <Dimension> <Steps> <Output Name> <Seed>"<<std::endl;
    return -1;
  }

  seed_rand(std::stoul(argv[4]));

  std::cout<<"A FastScape RB+PQ+GPU"<<std::endl;
  std::cout<<"C Richard Barnes TODO"<<std::endl;
  std::cout<<"h git_hash    = "<<GIT_HASH<<std::endl;
  std::cout<<"m Random seed = "<<argv[4]<<std::endl;

  const int width  = std::stoi(argv[1]);
  const int height = std::stoi(argv[1]);
  const int nstep  = std::stoi(argv[2]);

  CumulativeTimer tmr(true);
  const auto h = FastScape(width,height,nstep);
  std::cout<<"t Total calculation time    = "<<std::setw(15)<<tmr.elapsed()<<" microseconds"<<std::endl;

  PrintDEM(argv[3], h, width, height);

  delete[] h;

  return 0;
}
