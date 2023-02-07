#include <bits/stdc++.h>

#include "timer.h"
#include "int512.hpp"

using namespace std;

#define DONT_OPTIMIZE(x) { void* volatile dummy = &x; }

void fib_test() {
    Timer t("test_addition_fib");    
    for(int j = 0; j < 1000; j++){
        uint512_t a = 1, b = 1, c;
        DONT_OPTIMIZE(a);
        DONT_OPTIMIZE(b);
        DONT_OPTIMIZE(c);
        // fibonacci numbers
        for(int i = 0; i < 650; i++){
            c = a + b;
            a = b;
            b = c;
        }
    }
}

void gauss_test() {
    Timer t("test_addition_gauss");    
    uint512_t sum = 0;
    DONT_OPTIMIZE(sum);
    const int n = 100000000;
    for(int i = 1; i <= n; i++){
        sum = sum + i;
    }
    if (sum+sum != uint512_t(n)*uint512_t(n+1)){
        cout << "wrong result" << endl;
    }
}

void factorial_test() {
  Timer t("test_multiplication_factorial");    
  for(int j = 0; j < 50000; j++){
    uint512_t prod = 1;
    DONT_OPTIMIZE(prod);
    const int n = 90;
    for(int i = 1; i <= n; i++){
        prod = prod * i;
    }
  }
}

void collatz_test() {
  Timer t("test_division_collatz");    
  uint512_t n = "338838125384605298683130670138635876715138828490733673519118221753294017425139038491279";
  DONT_OPTIMIZE(n);
  uint512_t two = 2;
  while(n != 1){
    if (n % two == 0){
      n = n / two;
    } else {
      n = n * 3 + 1;
    }
  }
}

void mersenne_prime_test() {
  Timer t("test_division_mersenne_prime");
  //Find the first n mersenne prime
  const int n = 8;
  auto is_prime = [](const uint512_t& p) -> bool{
    if(p == 1) return false;
    if(p == 2) return true;
    if(p % 2 == 0) return false;
    uint512_t i = 3;
    while(i * i <= p){
      if(p - p / i * i == 0) return false;
      i = i + 2;
    }
    return true;
  };


  uint512_t p = 3;
  uint512_t two_to_i = 2;
  int twop = 1;
  DONT_OPTIMIZE(p);
  for(int i = 1; i <= n;){
    two_to_i = two_to_i * 2;
    twop++;
    p = two_to_i - 1;
    if(is_prime(p)){
      i++;
      cout << "Found prime " << i << ": 2^" << twop << "-1 = "  << p << endl;
    }

  }
}

auto factorize(long long n) -> vector<long long> {
  vector<long long> factors;
  for (long long i = 2; i * i <= n; i++) {
    while (n % i == 0) {
      factors.push_back(i);
      n /= i;
    }
  }
  if (n > 1)
    factors.push_back(n);
  return factors;
}


int main() {
 
  fib_test();
  gauss_test();
  factorial_test();
  collatz_test();
  mersenne_prime_test();
  return 0;
}