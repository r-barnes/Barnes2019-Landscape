CFLAGS = -O3 -march=native -g -fopt-info -fopt-info-vec-missed
WARNINGS = -Wall -Wpedantic -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef

.PHONY: all

all: fastscape_BW.exe fastscape_BW+P.exe fastscape_RB.exe fastscape_RB+P.exe fastscape_RB+PQ.exe

fastscape_BW.exe: fastscape_BW.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_BW.exe Timer.cpp CumulativeTimer.cpp fastscape_BW.cpp -Wno-unknown-pragmas   

fastscape_BW+P.exe: fastscape_BW+P.cpp
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_BW+P.exe Timer.cpp CumulativeTimer.cpp fastscape_BW+P.cpp -fopenmp

fastscape_RB.exe: fastscape_RB.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB.exe Timer.cpp CumulativeTimer.cpp fastscape_RB.cpp -Wno-unknown-pragmas   

fastscape_RB+P.exe: fastscape_RB+P.cpp  
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+P.exe Timer.cpp CumulativeTimer.cpp fastscape_RB+P.cpp -fopenmp

fastscape_RB+PQ.exe: fastscape_RB+PQ.cpp	
	$(CXX) $(CFLAGS) $(WARNINGS) -o fastscape_RB+PQ.exe Timer.cpp CumulativeTimer.cpp fastscape_RB+PQ.cpp	-fopenmp

clean:
	rm -rf *.exe