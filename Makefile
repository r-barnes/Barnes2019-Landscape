CFLAGS = -O3 -march=native -g -fopenmp
WARNINGS = -Wall -Wpedantic -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef

all: serial openmp

serial:
	g++ $(CFLAGS) $(WARNINGS) -o serial.exe fastscape_stack_serial2.cpp

openmp:
	g++ $(CFLAGS) $(WARNINGS) -o openmp.exe fastscape_openmp.cpp