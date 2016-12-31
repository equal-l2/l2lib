#pragma once
#include <chrono>

namespace l2lib{
  class time_counter{
    public:
      void start(){t1 = std::chrono::steady_clock::now();}
      void stop(){t2 = std::chrono::steady_clock::now();}
      template<typename Duration = std::chrono::milliseconds>
      typename Duration::rep count(){
        return std::chrono::duration_cast<Duration>(t2-t1).count();
      }
    private:
      std::chrono::time_point<std::chrono::steady_clock> t1;
      std::chrono::time_point<std::chrono::steady_clock> t2;
  };
}
