cmake_minimum_required (VERSION 3.9)

project(barnes2019_lem
  DESCRIPTION "High-performance landscape evolution models"
  LANGUAGES CXX CUDA
)

find_package(OpenMP REQUIRED)
find_package(OpenACC)

# GIT_HASH=`git rev-parse HEAD`
# COMPILE_TIME=`date -u +'%Y-%m-%d %H:%M:%S UTC'`

# CFLAGS = -O3 -march=native -g -DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" #-fopt-info -fopt-info-vec-missed  #-ftree-vectorize -funsafe-math-optimizations
# WARNINGS = -Wall -Wpedantic -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel -Wswitch-default -Wundef

# .PHONY: all

SET(COMMON_SOURCES
  CumulativeTimer.cpp
  random.cpp
)

SET(GIT_HASH "\"Unknown\"")

add_executable(fastscape_BW.exe
	fastscape_BW.cpp ${COMMON_SOURCES}
)
target_compile_options(fastscape_BW.exe PRIVATE -Wall -Wextra)
target_compile_definitions(fastscape_BW.exe PRIVATE GIT_HASH=${GIT_HASH})

add_executable(fastscape_BW+P.exe
	fastscape_BW+P.cpp ${COMMON_SOURCES}
)
target_compile_options(fastscape_BW+P.exe PRIVATE -Wall -Wextra)
target_compile_definitions(fastscape_BW+P.exe PRIVATE GIT_HASH=${GIT_HASH})
target_link_libraries(fastscape_BW+P.exe PRIVATE OpenMP::OpenMP_CXX)

add_executable(fastscape_BW+PI.exe
	fastscape_BW+PI.cpp ${COMMON_SOURCES}
)
target_compile_options(fastscape_BW+PI.exe PRIVATE -Wall -Wextra)
target_compile_definitions(fastscape_BW+PI.exe PRIVATE GIT_HASH=${GIT_HASH})
target_link_libraries(fastscape_BW+PI.exe PRIVATE OpenMP::OpenMP_CXX)

add_executable(fastscape_RB.exe
	fastscape_RB.cpp ${COMMON_SOURCES}
)
target_compile_options(fastscape_RB.exe PRIVATE -Wall -Wextra)
target_compile_definitions(fastscape_RB.exe PRIVATE GIT_HASH=${GIT_HASH})

add_executable(fastscape_RB+P.exe
	fastscape_RB+P.cpp ${COMMON_SOURCES}
)
target_compile_options(fastscape_RB+P.exe PRIVATE -Wall -Wextra)
target_compile_definitions(fastscape_RB+P.exe PRIVATE GIT_HASH=${GIT_HASH})
target_link_libraries(fastscape_RB+P.exe PRIVATE OpenMP::OpenMP_CXX)

add_executable(fastscape_RB+PI.exe
	fastscape_RB+PI.cpp ${COMMON_SOURCES}
)
target_compile_options(fastscape_RB+PI.exe PRIVATE -Wall -Wextra)
target_compile_definitions(fastscape_RB+PI.exe PRIVATE GIT_HASH=${GIT_HASH})
target_link_libraries(fastscape_RB+PI.exe PRIVATE OpenMP::OpenMP_CXX)

add_executable(fastscape_RB+PQ.exe
	fastscape_RB+PQ.cpp ${COMMON_SOURCES}
)
target_compile_options(fastscape_RB+PQ.exe PRIVATE -Wall -Wextra)
target_compile_definitions(fastscape_RB+PQ.exe PRIVATE GIT_HASH=${GIT_HASH})
target_link_libraries(fastscape_RB+PQ.exe PRIVATE OpenMP::OpenMP_CXX)

if(OpenACC_CXX_FOUND)
add_executable(fastscape_RB+GPU.exe
	fastscape_RB+GPU.cpp ${COMMON_SOURCES}
)
target_compile_options(fastscape_RB+GPU.exe PRIVATE -Wall -Wextra)
target_compile_definitions(fastscape_RB+GPU.exe PRIVATE GIT_HASH=${GIT_HASH})
target_link_libraries(fastscape_RB+GPU.exe PRIVATE OpenACC::OpenACC_CXX)
else()
  message(WARNING "OpenACC not found! Skipping fastscape_RB+GPU.exe")
endif()

	# echo "\033[91mCompiling 'fastscape_RB+GPU.exe' without OpenACC. No GPU acceleration will be used.\033[39m"
