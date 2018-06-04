//This file contains a serial implementation of Braun and Willett's FastScape
//algorithm. The implementation was developed by adapting Fortran code provided
//by Braun and attempts to be a faithful reproduction of the ideas therein.
#ifndef _quickscape_bw_hpp_
#define _quickscape_bw_hpp_

#include <array>
#include <vector>
#include "CumulativeTimer.hpp"

///The entire model is contained in a handy class, which makes it easy to set up
///and solve many such models.
class FastScape_BW {
 private:
  //Value used to indicate that a cell had no downhill neighbour and, thus, does
  //not flow anywhere.
  const int    NO_FLOW = -1;
  const double SQRT2   = 1.414213562373095048801688724209698078569671875376948; //Yup, this is overkill.
  const double dr[8]   = {1,SQRT2,1,SQRT2,1,SQRT2,1,SQRT2}; //Distance between adjacent cell centers on a rectangular grid arbitrarily scale to cell edge lengths of 1

 public:
  //NOTE: Having these constants specified in the class rather than globally
  //results in a significant speed loss. However, it is better to have them here
  //under the assumption that they'd be dynamic in a real implementation.
  double keq         = 2e-6;   //Stream power equation constant (coefficient)
  double neq         = 2;      //Stream power equation constant (slope modifier)
  double meq         = 0.8;    //Stream power equation constant (area modifier)
  double ueq         = 2e-3;   //Rate of uplift
  double dt          = 1000.;  //Timestep interval
  double tol         = 1e-3;   //Tolerance for Newton-Rhapson convergence while solving implicit Euler
  double cell_area   = 40000;  //Area of a single cell

 private:
  int width;        //Width of DEM
  int height;       //Height of DEM
  int size;         //Size of DEM (width*height)

  //The dataset is set up so that the outermost edge is never actually used:
  //it's a halo which allows every cell that is actually processed to consider
  //its neighbours in all 8 directions without having to check first to see if
  //it is an edge cell. The ring of second-most outer cells is set to a fixed
  //value to which everything erodes (in this model).

  //Rec directions (also used for nshift offsets) - see below for details
  //1 2 3
  //0   4
  //7 6 5

  std::vector<double> h;        //Digital elevation model (height)
  std::vector<double> accum;    //Flow accumulation at each point
  std::vector<int>    rec;      //Direction of receiving cell
  std::vector<int>    donor;    //Indices of a cell's donor cells
  std::vector<int>    ndon;     //How many donors a cell has
  std::vector<int>    stack;    //Indices of cells in the order they should be processed
  std::array<int,8>   nshift;   //Offset from a focal cell's index to its neighbours in terms of flat indexing

  //Timers for keeping track of how long each part of the code takes
  CumulativeTimer Tmr_Step1_Initialize;
  CumulativeTimer Tmr_Step2_DetermineReceivers;
  CumulativeTimer Tmr_Step3_DetermineDonors;
  CumulativeTimer Tmr_Step4_GenerateOrder;
  CumulativeTimer Tmr_Step5_FlowAcc;
  CumulativeTimer Tmr_Step6_Uplift;
  CumulativeTimer Tmr_Step7_Erosion;
  CumulativeTimer Tmr_Overall;

 public:
  ///Initializing code
  FastScape_BW(const int width0, const int height0);

 private:
  ///The receiver of a focal cell is the cell which receives the focal cells'
  ///flow. Here, we model the receiving cell as being the one connected to the
  ///focal cell by the steppest gradient. If there is no local gradient, than
  ///the special value NO_FLOW is assigned.
  void ComputeReceivers();

  ///The donors of a focal cell are the neighbours from which it receives flow.
  ///Here, we identify those neighbours by inverting the Receivers array.
  void ComputeDonors();

  void FindStack(const int c, int &nstack);

  ///Cells must be ordered so that they can be traversed such that higher cells
  ///are processed before their lower neighbouring cells. This method creates
  ///such an order.
  void GenerateOrder();

  ///Compute the flow accumulation for each cell: the number of cells whose flow
  ///ultimately passes through the focal cell multiplied by the area of each
  ///cell. Each cell could also have its own weighting based on, say, average
  ///rainfall.
  void ComputeFlowAcc();

  ///Raise each cell in the landscape by some amount, otherwise it wil get worn
  ///flat (in this model, with these settings)
  void AddUplift();

  ///Decrease he height of cells according to the stream power equation; that
  ///is, based on a constant K, flow accumulation A, the local slope between
  ///the cell and its receiving neighbour, and some judiciously-chosen constants
  ///m and n.
  ///    h_next = h_current - K*dt*(A^m)*(Slope)^n
  ///We solve this equation implicitly to preserve accuracy
  void Erode();

 public:

  ///Run the model forward for a specified number of timesteps. No new
  ///initialization is done. This allows the model to be stopped, the terrain
  ///altered, and the model continued. For space-efficiency, a number of
  ///temporary arrays are created each time this is run, so repeatedly running
  ///this function for the same model will likely not be performant due to
  ///reallocations. If that is your use case, you'll want to modify your code
  ///appropriately.
  void run(const int nstep);

  ///Returns a pointer to the data so that it can be copied, printed, &c.
  double* getH();
};

#endif
