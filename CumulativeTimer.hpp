#ifndef _timer_hpp_
#define _timer_hpp_

#include <chrono>
#include <cstdint>

class CumulativeTimer {
 private:
  typedef std::chrono::high_resolution_clock clock;
  typedef std::chrono::duration<double, std::ratio<1> > second;
  std::chrono::time_point<clock> start_time;
  uint64_t cumulative_time = 0;
  bool     running = false;
 public:
  CumulativeTimer  ();
  void     start   ();
  void     stop    ();
  void     reset   ();
  uint64_t elapsed () const;
};

#endif
