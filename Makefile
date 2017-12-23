GIT_HASH=`git rev-parse HEAD`
COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`

CFLAGS = -O3 -march=native -g -DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" #-fopt-info -fopt-info-vec-missed 
WARNINGS = -Wall -Wpedantic -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef

.PHONY: all

all: fastscape_BW.exe fastscape_BW+P.exe fastscape_RB.exe fastscape_RB+P.exe fastscape_RB+PQ.exe fastscape_RB+GPU.exe

fastscape_BW.exe: fastscape_BW.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_BW.exe    CumulativeTimer.cpp  random.cpp  fastscape_BW.cpp        -Wno-unknown-pragmas   

fastscape_BW+P.exe: fastscape_BW+P.cpp
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_BW+P.exe  CumulativeTimer.cpp  random.cpp  fastscape_BW+P.cpp      -fopenmp

fastscape_RB.exe: fastscape_RB.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB.exe    CumulativeTimer.cpp  random.cpp  fastscape_RB.cpp        -Wno-unknown-pragmas   

fastscape_RB+P.exe: fastscape_RB+P.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+P.exe  CumulativeTimer.cpp  random.cpp  fastscape_RB+P.cpp      -fopenmp

fastscape_RB+PQ.exe: fastscape_RB+PQ.cpp	
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+PQ.exe CumulativeTimer.cpp  random.cpp  fastscape_RB+PQ.cpp     -fopenmp

fastscape_RB+GPU.exe: fastscape_RB+GPU.cpp	
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+GPU.exe CumulativeTimer.cpp  random.cpp  fastscape_RB+GPU.cpp     -fopenmp

clean:
	rm -rf *.exe