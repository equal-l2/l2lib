#pragma once
#include <algorithm> // std::min
#include <cmath>     // std::pow
#include <iomanip>   // std::setw,setprecision
#include <iostream>  // std::clog
#include <sstream>   // std::ostringstream

namespace l2lib{
  class progress_bar{
    using big_uint = unsigned long long;
    const big_uint max;
    const unsigned width;
    const int precision;
    big_uint count;
    big_uint next_show_count;
    std::ostream& os;

    void show(){
      std::ostringstream oss;

      /* draw progress bar */
      oss << "\r[";

      const unsigned tics = static_cast<unsigned>((static_cast<double>(count)/max)*width);
      if(tics != 0){
        for(unsigned i = tics-1; i; --i){
          oss << '=';
        }
        oss << '>';
      }
      if(tics < width){
        for(unsigned i = width-tics; i; --i){
          oss << ' ';
        }
      }

      oss << "] [";

      /* draw progress in percentage */
      // `precision` > 0 : increased digits after decimal point
      // `precision` <= 0: print just before decimal point
      // if `precision` is negative, screen updating will happen less
      next_show_count = std::min({
          static_cast<big_uint>(static_cast<double>(tics+1)/width*max),                // when tics changes
          static_cast<big_uint>(count+max/(precision>0?std::pow(10,2+precision):100)), // when percentage changes
          max-1                                                                        // before the max
      });

      oss << std::fixed << std::setw(precision>0?4+precision:3) << std::setprecision(precision>0?precision:0)
          << count*(100.0/max) << "%]";

      // actual flushing of buffer will happen only when either progress in percentage or tics changes
      // `next_show_count` is capped with `max` to ensure last value is 100%

      os << oss.str() << std::flush;
    }
    public:
    progress_bar(big_uint max_in, unsigned width_in, int precision_in = 0, std::ostream& os_in = std::clog)
      // max       : max count (count value on 100% progress)
      // width     : width of progress bar in character
      // precision : digits after decimal point

      // Caution
      // Too much precision can cause great performance decrease

      :max(max_in),width(width_in),precision(precision_in),count(),next_show_count(),os(os_in){
        show();
      }
    ~progress_bar(){
      os << std::endl;
    }
    void adv(big_uint i = 1){
      if((count += i) >= next_show_count){show();}
    }
  };
}
