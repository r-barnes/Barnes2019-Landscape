#include "fastscape_BW.hpp"
#include "fastscape_BW+P.hpp"
#include "fastscape_BW+PI.hpp"
#include "fastscape_RB+GPU-openmp.hpp"
#include "fastscape_RB+GPU-graph.hpp"
#include "fastscape_RB.hpp"
#include "fastscape_RB+P.hpp"
#include "fastscape_RB+PI.hpp"
#include "fastscape_RB+PQ.hpp"
#include "random.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

///This is a quick-and-dirty, zero-dependency function for saving the outputs of
///the model in ArcGIS ASCII DEM format (aka Arc/Info ASCII Grid, AAIGrid).
///Production code for experimentation should probably use GeoTIFF or a similar
///format as it will have a smaller file size and, thus, save quicker.
void PrintDEM(
  const std::string filename, 
  const double *const h,
  const int width,
  const int height
){
  std::ofstream fout(filename.c_str());
  //Since the outer ring of the dataset is a halo used for simplifying
  //neighbour-finding logic, we do not save it to the output here.
  fout<<"ncols "<<(width- 2)<<"\n";
  fout<<"nrows "<<(height-2)<<"\n";
  fout<<"xllcorner 637500.000\n"; //Arbitrarily chosen value
  fout<<"yllcorner 206000.000\n"; //Arbitrarily chosen value
  fout<<"cellsize 500.000\n";     //Arbitrarily chosen value
  fout<<"NODATA_value -9999\n";   //Value which is guaranteed not to correspond to an actual data value
  for(int y=1;y<height-1;y++){
    for(int x=1;x<width-1;x++)
      fout<<h[y*width+x]<<" ";
    fout<<"\n";
  }
}



void GenerateRandomTerrain(double *const h, const int width, const int height){
  //srand(std::random_device()());
  for(int y=0;y<height;y++)
  for(int x=0;x<width;x++){
    const int c = y*width+x;
    h[c]  = uniform_rand_real(0,1);

    //Outer edge is set to 0 and never touched again. It is used only as a
    //convenience so we don't have to worry when a focal cell looks at its
    //neighbours.
    if(x == 0 || y==0 || x==width-1 || y==height-1)
      h[c] = 0;

    //Second outer-most edge is set to 0 and never touched again. This is the
    //baseline to which all cells would erode where it not for uplift. You can
    //think of this as being "sea level".
    if(x == 1 || y==1 || x==width-2 || y==height-2)
      h[c] = 0;
  }
}  



template<class T>
void RunModel(const int size, const int nstep, const std::string output_name, const std::string model_name){
  //Uses the RichDEM machine-readable line prefixes
  //Name of algorithm
  std::cout<<"A "<<model_name<<std::endl;           

  CumulativeTimer tmr(true);
  T tm(size,size);

  CumulativeTimer tmr_rand_gen(true);
  GenerateRandomTerrain(tm.getH(), size, size);
  tmr_rand_gen.stop();


  tm.run(nstep);
  std::cout<<"t Terrain generation time   = "<<std::setw(15)<<tmr_rand_gen.elapsed()<<" microseconds"<<std::endl;
  std::cout<<"t Total calculation time    = "<<std::setw(15)<<tmr.elapsed()         <<" microseconds"<<std::endl;

  PrintDEM(output_name, tm.getH(), size, size);
}



int main(int argc, char **argv){
  //Enable this to stop the program if a floating-point exception happens
  //feenableexcept(FE_ALL_EXCEPT);

  if(argc!=6){
    std::cerr<<"Syntax: "<<argv[0]<<" <Model> <Dimension> <Steps> <Output Name> <Seed>"<<std::endl;
    return -1;
  }

  const std::string model       = argv[1];
  const int         size        = std::stoi (argv[2]);
  const int         nstep       = std::stoi (argv[3]);
  const std::string output_name =            argv[4] ;
  const auto        rand_seed   = std::stoul(argv[5]);

  seed_rand(rand_seed);

  //Citation for algorithm
  std::cout<<"C Richard Barnes TODO"<<std::endl;
  //Git hash of code used to produce outputs of algorithm
  // std::cout<<"h git_hash    = "<<GIT_HASH<<std::endl;
  //Random seed used to produce outputs
  std::cout<<"m Random seed = "<<rand_seed<<std::endl;

  // if     (model=="BW")          RunModel<FastScape_BW>   (size, nstep, output_name, "fastscape_BW"    );
  // else if(model=="BW+P")        RunModel<FastScape_BWP>  (size, nstep, output_name, "fastscape_BW+P"  );
  // else if(model=="BW+PI")       RunModel<FastScape_BWPI> (size, nstep, output_name, "fastscape_BW+PI" );
  // else if(model=="RB")          RunModel<FastScape_RB>   (size, nstep, output_name, "fastscape_RB"    );
  // else if(model=="RB+P")        RunModel<FastScape_RBP>  (size, nstep, output_name, "fastscape_RB+P"  );
  // else if(model=="RB+PI")       RunModel<FastScape_RBPI> (size, nstep, output_name, "fastscape_RB+PI" );
  // else if(model=="RB+PQ")       RunModel<FastScape_RBPQ> (size, nstep, output_name, "fastscape_RB+PQ" );
  if(model=="RB+GPU")      RunModel<FastScape_RBGPU>(size, nstep, output_name, "fastscape_RB+GPU");
  // else if(model=="RB+GPUgraph") RunModel<FastScape_RBGPUgraph>(size, nstep, output_name, "fastscape_RB+GPUgraph");
  else {
    std::cerr<<"Unknown model! Choices: BW, BW+P, BW+PI, RB, RB+P, RB+PI, RB+PQ, RB+GPU, RB+GPUgraph."<<std::endl;
    return -1;
  }

  return 0;
}
