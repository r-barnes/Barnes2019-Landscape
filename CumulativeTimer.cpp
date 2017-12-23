#include "CumulativeTimer.hpp"
#include <stdexcept>

CumulativeTimer::CumulativeTimer(bool started){
  if(started)
    start();
}

void CumulativeTimer::start(){
  #pragma omp master
  {
    running    = true;
    start_time = clock::now();
  }
}

void CumulativeTimer::stop(){
  #pragma omp master
  {
    if(!running)
      throw std::runtime_error("Can't stop a Timer that hasn't been started!");

    running          = false;
    cumulative_time += static_cast<unsigned long int>(std::chrono::duration_cast<std::chrono::microseconds> (clock::now() - start_time).count());
  }
}

void CumulativeTimer::reset(){
  cumulative_time = 0;
}

uint64_t CumulativeTimer::elapsed() const { 
  if(running)
    return cumulative_time + static_cast<unsigned long int>(std::chrono::duration_cast<std::chrono::microseconds> (clock::now() - start_time).count()); 
  else
    return cumulative_time;
}
