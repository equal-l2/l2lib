#pragma once
#include <complex>
#include <vector>

namespace l2lib{
  constexpr long double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286L;
  using real_t = long double;
  using comp_t = std::complex<real_t>;
  using data_vector = std::vector<comp_t>;
  data_vector dft(const data_vector& v){
    using namespace std::literals::complex_literals;
    if(v.size() == 0) return {};
    const auto n = v.size();
    const real_t fac = 2*pi/n;
    data_vector ret(n,0);
    for(size_t i=0; i < n; ++i){
      comp_t val = 0;
      for(size_t j=0; j < n; ++j){
        val += v[j]*std::exp(-1il*(fac*i*j));
      }
      ret[i] = val;
    }
    return ret;
  }

  data_vector idft(data_vector v){
    if(v.size() == 0) return {};
    const auto n = v.size();
    for(auto &val : v){
      val = conj(val);
    }
    v = dft(v);
    for(auto &val : v){
      val = conj(val)/real_t(n);
    }
    return v;
  }

  #if L2LIB_DFT_FFTW_ENABLE
  data_vector dft_fftw(const data_vector& v){
    if(v.size() == 0) return {};
    const auto n = v.size();
    data_vector ret(n,0);
    auto plan = fftwl_plan_dft_1d(n,(fftwl_complex*)v.data(),(fftwl_complex*)ret.data(),FFTW_FORWARD,FFTW_ESTIMATE);
    fftwl_execute(plan);
    return ret;
  }

  data_vector idft_fftw(data_vector& v){
    if(v.size() == 0) return {};
    const auto n = v.size();
    data_vector ret(n,0);
    auto plan = fftwl_plan_dft_1d(n,(fftwl_complex*)v.data(),(fftwl_complex*)ret.data(),FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwl_execute(plan);
    for(auto& v : ret){
      v /= n;
    }
    return ret;
  }
  #endif
}
