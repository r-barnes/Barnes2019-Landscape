#ifndef _timer_hpp_
#define _timer_hpp_

#include <chrono>

class Timer {
 private:
  typedef std::chrono::high_resolution_clock clock;
  typedef std::chrono::duration<double, std::ratio<1> > second;
  std::chrono::time_point<clock> start_time;
  std::chrono::time_point<clock> old_time;
 public:
  Timer();
  void   reset  ();
  double elapsed() const;
  double lap    ();
};

#endif
