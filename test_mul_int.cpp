#include "mul_int.hpp"
#include "time_counter.hpp"
#include <iostream>
#include <cstdlib>
#include <random>
#include <functional>

#define OP_TEST_BASE(op, ex_zero)\
{\
  std::cout << "start test for operator " << #op << '\n';\
  l2lib::time_counter tc;\
  tc.start();\
  for(int i=0; i<1000000;++i){\
    long long l1,l2;\
    l2lib::mul_int m1,m2;\
    l1 = gen();\
    l2 = gen();\
    if(ex_zero && l2 == 0) continue;\
    m1 = l1;\
    m2 = l2;\
    if(l1 op l2 != m1 op m2){\
      std::cout << "failed test for operator " << #op << '\n';\
      std::cout << "values  : " << l1 << " " << l2 << '\n';\
      std::cout << "expected: " << (l1 op l2) << '\n';\
      std::cout << "actual  : " << (m1 op m2) << '\n';\
      exit(1);\
    }\
  }\
  tc.stop();\
  std::cout << tc.count<>() << std::endl;\
}

#define OP_TEST(op) OP_TEST_BASE(op,false)
#define OP_TEST_EX_ZERO(op) OP_TEST_BASE(op,true)

int main(){
  auto gen = std::bind(std::uniform_int_distribution<long long>(-10000,10000),std::mt19937());
  OP_TEST(+);
  OP_TEST(-);
  OP_TEST(*);
  OP_TEST_EX_ZERO(/);
  OP_TEST_EX_ZERO(%);
}
