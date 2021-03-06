GIT_HASH=`git rev-parse HEAD`
COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`

CFLAGS = -O3 -march=native -g -DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" #-fopt-info -fopt-info-vec-missed  #-ftree-vectorize -funsafe-math-optimizations
WARNINGS = -Wall -Wpedantic -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef

.PHONY: all

all: fastscape_BW.exe fastscape_BW+P.exe fastscape_BW+PI.exe fastscape_RB.exe fastscape_RB+P.exe fastscape_RB+PI.exe fastscape_RB+PQ.exe fastscape_RB+GPU.exe

fastscape_BW.exe: fastscape_BW.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_BW.exe    CumulativeTimer.cpp  random.cpp  fastscape_BW.cpp        -Wno-unknown-pragmas   

fastscape_BW+P.exe: fastscape_BW+P.cpp
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_BW+P.exe  CumulativeTimer.cpp  random.cpp  fastscape_BW+P.cpp      -fopenmp

fastscape_BW+PI.exe: fastscape_BW+PI.cpp
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_BW+PI.exe  CumulativeTimer.cpp  random.cpp  fastscape_BW+PI.cpp      -fopenmp

fastscape_RB.exe: fastscape_RB.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB.exe    CumulativeTimer.cpp  random.cpp  fastscape_RB.cpp        -Wno-unknown-pragmas   

fastscape_RB+P.exe: fastscape_RB+P.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+P.exe  CumulativeTimer.cpp  random.cpp  fastscape_RB+P.cpp      -fopenmp

fastscape_RB+PI.exe: fastscape_RB+PI.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+PI.exe  CumulativeTimer.cpp  random.cpp  fastscape_RB+PI.cpp      -fopenmp

fastscape_RB+PQ.exe: fastscape_RB+PQ.cpp	
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+PQ.exe CumulativeTimer.cpp  random.cpp  fastscape_RB+PQ.cpp     -fopenmp

fastscape_RB+GPU.exe: fastscape_RB+GPU.cpp
	echo "\033[91mCompiling 'fastscape_RB+GPU.exe' without OpenACC. No GPU acceleration will be used.\033[39m"
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+GPU.exe CumulativeTimer.cpp  random.cpp  fastscape_RB+GPU.cpp -Wno-unknown-pragmas -Wno-shadow

clean:
	rm -rf *.exe