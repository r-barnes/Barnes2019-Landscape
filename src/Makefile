GIT_HASH=`git rev-parse HEAD`
COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`

CFLAGS = -O3 -march=native -g -DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" #-fopt-info -fopt-info-vec-missed  #-ftree-vectorize -funsafe-math-optimizations
WARNINGS = -Wall -Wpedantic -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef

.PHONY: all

all: quickscape.exe

quickscape.exe: random.cpp CumulativeTimer.cpp fastscape_BW.cpp fastscape_BW+P.cpp fastscape_BW+PI.cpp fastscape_RB.cpp fastscape_RB+P.cpp fastscape_RB+PI.cpp fastscape_RB+PQ.cpp fastscape_RB+GPU.cpp fastscape_BW.hpp fastscape_BW+P.hpp fastscape_BW+PI.hpp fastscape_RB.hpp fastscape_RB+P.hpp fastscape_RB+PI.hpp fastscape_RB+PQ.hpp fastscape_RB+GPU.hpp random.hpp CumulativeTimer.hpp main.cpp
	echo "\033[91mCompiling 'fastscape_RB+GPU.exe' without OpenACC. No GPU acceleration will be used.\033[39m"
	$(CXX) $(CFLAGS) $(WARNINGS) -o quickscape.exe  main.cpp CumulativeTimer.cpp  random.cpp  fastscape_BW.cpp fastscape_BW+P.cpp fastscape_BW+PI.cpp fastscape_RB.cpp fastscape_RB+P.cpp fastscape_RB+PI.cpp fastscape_RB+PQ.cpp fastscape_RB+GPU.cpp -Wno-unknown-pragmas -fopenmp 

clean:
	rm -rf *.exe