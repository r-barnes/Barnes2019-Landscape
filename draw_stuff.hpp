// void DrawColours(){
//   for i in range(rmin,rmax+1):
//     print("\\definecolor{{color{v}}}{{RGB}}{{{gs},{gs},{gs}}}".format(v=i, gs=int(((i-rmin)*255)/rrange)))  
// }

  // StartPicture();
  // DrawStack(h.data(), WIDTH, HEIGHT, nshift.data(), rec.data(), stack.data());
  // EndPicture();

const int COLORNUM=12;

#include <map>
#include <algorithm>
#include <numeric>
#include <random>

const unsigned int RUN_NUMBER = std::random_device()();

void DrawDEM(
  const double *const elev,
  const int width,
  const int height,
  const bool fill
){
  if(width>=100)
    return;
  for(int y=0;y<height;y++)
  for(int x=0;x<width;x++){
    if(fill)
      std::cout<<"\\filldraw[draw=black] ("<<x<<","<<y<<") rectangle ("<<(x+1)<<","<<(y+1)<<");\n"; //,fill=color{v}
    else
      std::cout<<"\\draw ("<<x<<","<<y<<") rectangle ("<<(x+1)<<","<<(y+1)<<");\n";
  }
}

void DrawStack(
  const double *const elev,
  const int width,
  const int height,  
  const int    *const nshift,
  const int    *const rec,
  const int    *const stack
){
  int color     = -1;
  int prevstart = -1;

  DrawDEM(elev,width,height,false);

  std::vector<int> stacksizes;

  for(int i=0;i<width*height;i++){
    const int c     = stack[i];
    const double cy = (c/width)+0.5;
    const double cx = (c%width)+0.5;
    if(rec[c]==NO_FLOW){
      stacksizes.push_back(i-prevstart);
      prevstart = i;
      color++;
      continue;
    }    
    const int n     = c+nshift[rec[c]];
    const double ny = (n/width)+0.5;
    const double nx = (n%width)+0.5;
    // std::cout<<"\\draw[<-,line width=4, color="<<color<<"] ("<<cx<<","<<cy<<") -- ("<<nx<<","<<ny<<");\n";
    // std::cout<<"\\draw[fill="<<color<<"] ("<<cx<<","<<cy<<") circle [radius=0.25];\n";
    // std::cout<<"\\draw[fill="<<color<<"] ("<<nx<<","<<ny<<") circle [radius=0.25];\n";
    if(width<100){
      std::cout<<"\\draw[-{Latex[length=2.5mm, width=2mm]},line width=2,color=color"<<(color%COLORNUM)<<"] ("<<cx<<","<<cy<<") -- ("<<nx<<","<<ny<<");\n";
      std::cout<<"\\path[fill=color"<<(color%COLORNUM)<<"] ("<<cx<<","<<cy<<") circle [radius=0.125];\n";
      std::cout<<"\\path[fill=color"<<(color%COLORNUM)<<"] ("<<nx<<","<<ny<<") circle [radius=0.125];\n";  
    }  
  }

  std::map<int, int> ssmap;
  for(const auto x: stacksizes)
    ssmap[x]++;

  std::vector<int> vals;
  std::vector<int> count;
  for(auto kv=ssmap.rbegin();kv!=ssmap.rend();kv++){
    vals.push_back(kv->first);
    count.push_back(kv->second);
  }

  std::vector<int> cumsum(count.size());
  std::partial_sum(count.begin(), count.end(), cumsum.begin());

  std::reverse(vals.begin(), vals.end());
  std::reverse(cumsum.begin(), cumsum.end());
  for(int i=0;i<vals.size();i++)
    std::cout<<"stacksize "<<RUN_NUMBER<<" "<<vals.at(i)<<" "<<cumsum.at(i)<<std::endl;
}

void DrawQueue(
  const double *const elev,
  const int width,
  const int height,
  const int    *const nshift,
  const int    *const rec,
  const int    *const queue,
  const int    *const levels,
  const int nlevels
){
  DrawDEM(elev,width,height,false);

  std::vector<std::string> lnames = {{"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F"}};

  std::cerr<<"nlevels = "<<nlevels<<std::endl;

  for(int li=0;li<nlevels-1;li++){
    const int color = li;
    for(int i=levels[li];i<levels[li+1];i++){
      const int c     = queue[i];
      const double cy = c/width+0.5;
      const double cx = c%width+0.5;
      if(rec[c]==NO_FLOW)
        continue;
      const int n     = c+nshift[rec[c]];
      const double ny = n/width+0.5;
      const double nx = n%width+0.5;
      if(width<100)
        std::cout<<"\\draw[-{Latex[length=2.5mm, width=2mm]},line width=2] ("<<cx<<","<<cy<<") -- ("<<nx<<","<<ny<<");\n";
    }
  }
  for(int li=0;li<nlevels-1;li++){
    const int color = li;
    std::cout<<"levelsize "<<RUN_NUMBER<<" "<<li<<" "<<(levels[li+1]-levels[li])<<std::endl;
    for(int i=levels[li];i<levels[li+1];i++){
      const int c     = queue[i];
      const double cy = c/width+0.5;
      const double cx = c%width+0.5;
      if(rec[c]==NO_FLOW)
        continue;
      const int n     = c+nshift[rec[c]];
      const double ny = n/width+0.5;
      const double nx = n%width+0.5;
      if(width<100)
        std::cout<<"\\draw[fill=color"<<(li%COLORNUM)<<"] ("<<cx<<","<<cy<<") circle [radius=0.3] node {\\scriptsize "<<lnames[li]<<"};\n";
    }
  }  
}

void StartPicture(){
  std::cout<<"\\begin{tikzpicture}[scale=0.35]\n";
}

void EndPicture(){
  std::cout<<"\\end{tikzpicture}\n";
}