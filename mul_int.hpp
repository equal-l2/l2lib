#pragma once
#include <string>
#include <functional>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <sstream>

//#include "dft.hpp"

#include <iostream>
#include <iomanip>
#include <cassert>

#if BOUND_CHECK
  #define ACC(n) digits.at(n)
#else
  #define ACC(n) digits[n]
#endif

namespace l2lib{
  class mul_int{ // multiprecision integer
  public:
    using self = mul_int;
    using dig_t = uint_fast8_t;
    using dig_vec = std::vector<dig_t>;
    using num_t = long long;

  private:
    dig_vec digits;
    bool sgn; // true when the number is negative

  public:
    mul_int(const std::string s = "0"):digits(),sgn(false){
      if(s.empty()) return;

      const auto begin = s.rbegin();
      auto end = s.rend();
      if(s[0] == '-'){ // if s represents negative number
        sgn = true;
        --end;
      }

      while(*(end - 1) == '0' && end != begin) {
        // inline zero suppression
        --end;
      }

      if(end == begin) { // if s == "0" or "-0"
        sgn = false;
        digits.push_back(0);
        return;
      }

      digits.reserve(end-begin);
      for(auto it = begin; it != end; ++it){
        dig_t buf = ctoi(*it);
        digits.push_back(buf);
      }
    }

    mul_int(num_t i):digits(),sgn(i<0) {
      if(is_minus()) i = -i;

      do{
        digits.push_back(i%10);
        i /= 10;
      }while(i != 0);
    }

    mul_int(const self& n):digits(n.digits),sgn(n.sgn){}
    mul_int(const dig_vec& vec):digits(vec),sgn(false){}

    inline size_t size() const {
      // return size of number (number of digits)
      // not always same as 'real' size (size of vector digits)
      return zero_suppress()+1;
    }

    self& operator=(const self& rhs){
      digits = rhs.digits;
      sgn = rhs.sgn;
      return *this;
    }

    self& operator=(const num_t rhs){
      return operator=(self(rhs));
    }

    self& operator=(const std::string& rhs){
      return operator=(self(rhs));
    }

    inline self operator+() const {return *this;}

    self operator-() const {
      auto ret = *this;
      ret.flip_sgn();
      return ret;
    }

    self& operator+=(const self& rhs){
      if(rhs.is_minus()) return operator-=(-rhs);
      if(rhs.is_zero()) return *this;

      // guaranteed rhs > 0

      if(is_zero()){
        *this = rhs;
      }
      else{
        if(is_minus()){ // *this<0, rhs>0
          if(abs_lessthan(rhs)){ // |*this| < |rhs|
            auto rhs_copy = rhs;
            rhs_copy.abs_sub(*this);
            *this = rhs_copy;
          }
          else if(abs_equal(rhs)){ // |*this| == |rhs|
            *this = 0;
          }
          else{ // |*this| > |rhs|
            abs_sub(rhs);
          }
        }
        else{ // *this>0, rhs>0
          abs_add(rhs);
        }
      }

      return *this;
    }

    self& operator-=(const self& rhs){
      if(rhs.is_minus()) return operator+=(-rhs);
      if(rhs.is_zero()) return *this;

      // guaranteed rhs > 0

      if(is_zero()){
        *this = -rhs;
      }
      else{
        if(is_minus()){ // *this<0, rhs>0
          abs_add(rhs);
        }
        else{ // *this>0, rhs>0
          if(abs_lessthan(rhs)){ // |*this| < |rhs| -> *this-rhs < 0
            auto rhs_copy = rhs;
            rhs_copy.abs_sub(*this);
            *this = -rhs_copy;
          }
          else if(abs_equal(rhs)){ // |*this| == |rhs|
            *this = 0;
          }
          else{ // |*this| > |rhs|
            abs_sub(rhs);
          }
        }
      }

      return *this;
    }

    self& operator*=(const self& rhs){
      auto sgn_dif = (is_minus() != rhs.is_minus());
      if(is_zero() || rhs.is_zero()) *this = 0;
      else if(abs_is_ten_multiple()){
        auto shift_num = zero_suppress();
        *this = rhs;
        shift_up_digits(shift_num);
      }
      else if(rhs.abs_is_ten_multiple()) shift_up_digits(rhs.zero_suppress());
      else karatsuba_multiply(rhs);

      set_sgn(sgn_dif);

      return *this;
    }

    self& operator/=(const self& rhs){
      auto sgn_dif = (is_minus() != rhs.is_minus());
      if(rhs.is_zero()) throw std::runtime_error("zero division");
      if(is_zero() || abs_lessthan(rhs)) *this = 0;
      else if(abs_equal(rhs)) *this = 1;
      // guaranteed  |*this| > |rhs| != 0
      else if(!rhs.abs_is_one()){
        self cnt("");
        cnt.digits.resize(size(),0);
        auto offset = zero_suppress();
        for(auto i=size();i>0;--i){
          if(!abs_lessthan(rhs,offset)){
            cnt.ACC(i-1) = partial_div(rhs,offset);
          }
          --offset;
        }
        *this = cnt;
      }

      set_sgn(sgn_dif);
      return *this;
    }

    self& operator%=(const self& rhs){
      if(rhs.is_zero()) throw std::runtime_error("zero division");
      if(
          is_zero()
       || rhs.equal_to_one_dig(1)
       || (is_even() && rhs.equal_to_one_dig(2))
       || ((digits[0] == 0 || digits[0] == 5) && rhs.equal_to_one_dig(5))
       || abs_equal(rhs)
      ) *this = 0;
      // guaranteed |rhs| != (0,1,2,5), *this != 0, |*this| != |rhs|, *this%2 != 0, *this%5 != 0
      else if(!abs_lessthan(rhs)){ // *this > rhs
        auto offset = zero_suppress();
        for(auto i=size();i>0;--i){
          if(!abs_lessthan(rhs,offset)) partial_div(rhs,offset);
          --offset;
        }
        if(is_zero() && is_minus()) sgn = false;
      }
      return *this;
    }

    self& operator++(){
      return (*this += 1);
    }

    self operator++(int){
      auto ret = *this;
      operator++();
      return ret;
    }

    self& operator--(){
      return (*this -= 1);
    }

    self operator--(int){
      auto ret = *this;
      operator--();
      return ret;
    }

    friend std::ostream& operator<<(std::ostream& os, const self& val){
      if(val.digits.empty()) return (os << "(EMPTY)");
      return os << val.to_string();
    }

    friend std::istream& operator>>(std::istream& is, self& val){
      std::string buf;
      is >> buf;
      val = self(buf);
      return is;
    }

    #define AR_OP(op) \
    friend self operator op(self lhs, const self& rhs){\
      return (lhs op##= rhs);\
    }\

    AR_OP(+)
    AR_OP(-)
    AR_OP(*)
    AR_OP(/)
    AR_OP(%)

    #undef AR_OP

    #define REL_OP(op,ex) \
    friend bool operator op(const self& lhs, const self& rhs){\
      return ex; \
    } \
    \

    // TODO: optimize for 0
    REL_OP(==,(lhs.is_minus() == rhs.is_minus()) && lhs.abs_equal(rhs))
    REL_OP(<,(lhs.is_minus() && !rhs.is_minus()) || ((lhs.is_minus() == rhs.is_minus()) && lhs.abs_lessthan(rhs)))
    REL_OP(>,(rhs<lhs))
    REL_OP(!=,!(lhs == rhs))
    REL_OP(<=,!(lhs > rhs))
    REL_OP(>=,!(lhs < rhs))

    #undef REL_OP

    std::string to_string() const {
      std::ostringstream oss;
      if(is_minus()) oss << '-';
      auto top = size();
      for(auto i = top; i > 0; --i){
        oss << static_cast<unsigned>(get_digit(i-1));
      }
      return oss.str();
    }

    inline dig_t get_digit(size_t i) const { // zero-oriented
      return ((real_size() <= i) ? 0 : digits[i]);
    }

    /*
    self dft_mul(const self& rhs) const {
      const auto max = size()+rhs.size();

      std::vector<std::complex<long double>> a_dat;
      std::vector<std::complex<long double>> b_dat;

      for(size_t i=0; i<max; ++i){
        a_dat.push_back(get_digit(i));
        b_dat.push_back(rhs.get_digit(i));
      }

      a_dat = dft(a_dat);
      b_dat = dft(b_dat);

      std::transform(a_dat.begin(),a_dat.end(),b_dat.begin(),a_dat.begin(),std::multiplies<std::complex<long double>>());

      a_dat = idft(a_dat);

      self ret("");
      unsigned buf,carry = 0; // buf may contain larger value than max of dig_t
      for(const auto& v : a_dat){
        buf = carry+static_cast<unsigned>(v.real()+0.5L); // '0.5L' は四捨五入のため
        carry = buf/10;
        // guaranteed buf < 10  <=> containable with 8bit
        ret.digits.push_back(static_cast<dig_t>(buf%10));
      }
      if(!ret.is_zero() && (is_minus() != rhs.is_minus())) ret.flip_sgn();
      return ret;
    }
    */

  private:
    inline bool is_zero() const { return equal_to_one_dig(0); }
    inline bool abs_is_one() const { return equal_to_one_dig(1); }
    inline bool abs_is_ten_multiple() const { // check if number is in form of 100...(any number of 0)...000
      return get_digit(zero_suppress())==1 && 
          std::all_of(digits.begin(),digits.begin()+zero_suppress(),[](auto i){return i==0;});
    }

    // return if *this equals to n
    inline bool equal_to_one_dig(dig_t n) const { return (size() == 1) && (digits[0] == n); }

    inline bool is_minus() const { return sgn; }

    inline bool is_even() const {return (digits[0]%2 == 0);}
    inline size_t real_size() const {return digits.size();}

    inline void flip_sgn(){set_sgn(!sgn);}
    inline void set_sgn(bool sgn_arg){if(!is_zero()) sgn = sgn_arg;}

    inline std::size_t zero_suppress() const {
      auto top = real_size()-1;
      if(top != 0){ // if digits of number != 0
        for(;get_digit(top) == 0 && top > 0; --top);
      }
      return top;
    }

    self& long_multiply(const self& rhs){
      self buf;
      for(size_t i=0; i<size(); ++i){
        for(size_t j=0; j<rhs.size(); ++j){
          buf.digit_add(get_digit(i)*rhs.get_digit(j),i+j);
        }
      }
      *this = buf;
      return *this;
    }

    void karatsuba_multiply(const self& rhs){
      constexpr size_t div_num = 20;
      if(size() <= div_num && rhs.size() <= div_num) long_multiply(rhs);
      else{
        auto this_div = karatsuba_divide(div_num);
        auto rhs_div = rhs.karatsuba_divide(div_num);
        auto z0 = this_div.first*rhs_div.first;
        auto z2 = this_div.second*rhs_div.second;
        auto z1 = (this_div.first+this_div.second)*(rhs_div.first+rhs_div.second)-z0-z2;
        *this = z0 + z1.shift_up_digits(div_num) + z2.shift_up_digits(2*div_num);
      }
    }

    std::pair<self,self> karatsuba_divide(size_t div) const {
      return {slice(0,div-1),slice(div,zero_suppress())};
    }

    self& shift_up_digits(const size_t n){
      digits.insert(digits.begin(),n,0);
      return *this;
    }

    self slice(size_t begin_digit, size_t end_digit) const {
      if(begin_digit > end_digit) return self();
      if(begin_digit == end_digit) return self(dig_vec(1,get_digit(begin_digit)));
      dig_vec slice_vec(end_digit-begin_digit+1);
      for(size_t i=0;i<=end_digit-begin_digit;++i){
        slice_vec.at(i) = get_digit(begin_digit+i);
      }
      return self(slice_vec);
    }

    void check_carry(){
      for(size_t i=0; i<real_size(); ++i) digit_add(0,i);
    }

    void digit_add(const dig_t val,const size_t index) {
      //指定された桁に対して、valを足し算。繰上り考慮
      if(real_size()-1 < index) digits.push_back(0);
      const dig_t buf = get_digit(index) + val;
      ACC(index) = buf % 10;
      if(buf >= 10) digit_add(buf / 10,index + 1);
    }

    void digit_sub(const dig_t val,const size_t index) {
      //指定された桁に対して、valを引き算。繰下り考慮
      if(digits[index] < val) {
        ACC(index) += 10;
        digit_sub(1,index + 1);
      }
      ACC(index) -= val;
    }

    void abs_add(const self& rhs, const size_t offset = 0){
      for(size_t i=0;i<rhs.size();++i) digit_add(rhs.get_digit(i),i+offset);
    }

    void abs_sub(const self& rhs, const size_t offset = 0){
      // assume *this > rhs
      for(size_t i = 0;i < rhs.size();++i) digit_sub(rhs.get_digit(i),i+offset);
    }

    bool abs_equal(const self& rhs, const size_t offset = 0) const {
      if(size() < offset) return rhs.is_zero();
      auto this_size = size()-offset;
      if(this_size != rhs.size()) return false;
      for(size_t i=0;i<rhs.size();++i){
        if(get_digit(i+offset) != rhs.get_digit(i)) return false;
      }
      return true;
    }

    bool abs_lessthan(const self& rhs, const size_t offset = 0) const { // return true if *this < rhs
      if(size() < offset) return true;
      // guaranteed offset <= size()
      auto this_size = size() - offset;
      auto rhs_size = rhs.size();
      if(this_size < rhs_size) return true;
      if(this_size > rhs_size) return false;
      for(size_t i = rhs_size + 1;i > 0;--i) { // 上から
        if(get_digit(i - 1 + offset) < rhs.get_digit(i - 1)) return true;
        if(get_digit(i - 1 + offset) > rhs.get_digit(i - 1)) return false;
      }
      return false;
    }

    dig_t partial_div(const self& rhs, const size_t offset){
      dig_t cnt = 0;
      while(!abs_lessthan(rhs,offset)){ // while |part| >= |rhs|
        abs_sub(rhs,offset);
        ++cnt;
      }
      return cnt;
    }

    void hard_zero_suppress(){
      if(digits.empty()) return;
      const auto head = zero_suppress();
      if(head == real_size()-1) return;
      digits.erase(digits.begin()+1+static_cast<long long>(head),digits.end()-1);
    }

    inline dig_t ctoi(const char c) const {
      if('0' <= c && c <= '9') return static_cast<dig_t>(c-'0');
      else throw std::invalid_argument("'"+(std::string(1,c)+"' is out of range [0-9]"));
    }
  };

  mul_int abs(const mul_int& n){
    return (n < 0 ? -n : n);
  }

  mul_int pow(const mul_int& base, const mul_int& exp){
    if(exp == 0) return 1;
    if(exp == 1) return base;

    auto ret = base;
    for(mul_int i=1;i<exp;++i){
      ret *= base;
    }
    return ret;
  }

  mul_int factorial(const unsigned n){
    mul_int ret = 1;
    for(unsigned i=2; i<=n; ++i){
      ret *= i;
    }
    return ret;
  }

  bool is_prime(const mul_int& n){
    if(n <  2) return false;
    if(n == 2) return true;
    if(n == 3) return true;
    // guaranteed n > 3

    if(n%2 == 0) return false;
    if(n%3 == 0) return false;
    // guaranteed n is not divisible with 6

    const auto res = n%6;
    if(res != 1 && res != 5) return false;

    // 6の倍数前後の数を使って試し割りをする
    for(mul_int i=5;i*i<=n;i+=6){
      if(n%i     == 0) return false;
      if(n%(i+2) == 0) return false;
    }

    return true;
  }

  std::vector<mul_int> find_nontrivial_divisors(const mul_int& num){
    std::vector<mul_int> divisors;
    if(num%2 == 0){
      for(mul_int i=2; i<num; i++){
        if(num%i == 0) divisors.push_back(i);
      }
    }
    else{
      for(mul_int i=3; i<num; i+=2){
        if(num%i == 0) divisors.push_back(i);
      }
    }
    return divisors;
  }
}
