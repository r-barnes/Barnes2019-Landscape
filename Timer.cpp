#include "Timer.hpp"

Timer::Timer(){
  reset();
}

void Timer::reset(){
  start_time = clock::now();
  old_time   = start_time;
}

double Timer::elapsed() const { 
  return std::chrono::duration_cast<second> (clock::now() - start_time).count(); 
}

double Timer::lap() {
  const auto now  = clock::now();
  const auto diff = std::chrono::duration_cast<second> (now - old_time).count();
  old_time        = now;
  return diff;
}
