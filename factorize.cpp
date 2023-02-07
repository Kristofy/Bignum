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
        prod = prod * uint512_t(i);
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
      n = n * 3 + uint512_t(1);
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
  // __m256i MaskGeneratorLo = _mm256_set_epi32(9, 10, 11, 12, 13, 14, 15, 16);
  // __m256i MaskGeneratorHi = _mm256_set_epi32(1, 2, 3, 4, 5, 6, 7, 8);
    

  // for(int i = 0; i < 16; i++){
  //   __m256i MaskLo = _mm256_cmpgt_epi32(MaskGeneratorLo, _mm256_set1_epi32(i));
  //   __m256i MaskHi = _mm256_cmpgt_epi32(MaskGeneratorHi, _mm256_set1_epi32(i));
  //   __m256i A = _mm256_and_si256(MaskGeneratorLo, MaskLo);
  //   __m256i B = _mm256_and_si256(MaskGeneratorHi, MaskHi);
  //   int arr[16];
  //   _mm256_storeu_si256((__m256i *)arr, A);
  //   _mm256_storeu_si256((__m256i *)(arr+8), B);
  //   for(int j = 0; j < 16; j++){
  //     cout << arr[j] << " ";
  //   }
  //   cout << endl << endl;
  // }
  // return 0;

  // alignas(32) int a[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  // alignas(32) int b[8] = {9, 10, 11, 12, 13, 14, 15, 16};
  // static __m256i ShiftRMask = _mm256_set_epi32(6,5,4,3,2,1,0,7);

  // __m256i Zero = _mm256_setzero_si256();

  // __m256i X= _mm256_load_si256((__m256i *)a);
  // __m256i Y = _mm256_load_si256((__m256i *)b);

  // __m256i CL, CH;
  
  // Y = _mm256_insert_epi32(Y, _mm256_extract_epi32(X, 7), 7);
  // X = _mm256_insert_epi32(X, 0, 7);
  // CL = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(X), ShiftRMask));
  // CH = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(Y), ShiftRMask));


  // // store
  // _mm256_store_si256((__m256i *)a, CL);
  // _mm256_store_si256((__m256i *)b, CH);

  // // cout
  // for (int i = 0; i < 8; i++) {
  //   cout << a[i] << " ";
  // }
  // cout << endl;
  // for (int i = 0; i < 8; i++) {
  //   cout << b[i] << " ";
  // }
  // cout << endl;


  fib_test();
  gauss_test();
  factorial_test();
  collatz_test();
  mersenne_prime_test();
  return 0;
  vector<pair<long long, vector<long long>>> test_data = {{1223456789123LL, {}},   {1223456789126LL, {}},
                                                    {642524234131232LL, {}}, {1000000007LL, {}},
                                                    {49290915216859LL, {}},  {84272973607266299, {}}};
  chrono::steady_clock::time_point begin = chrono::steady_clock::now();
  for (auto &[n, factors] : test_data) {
    factors = factorize(n);
  }
  chrono::steady_clock::time_point end = chrono::steady_clock::now();

  for (auto &[n, factors] : test_data) {
    cout << n << ": ";
    for (auto &factor : factors) {
      cout << factor << " ";
    }
    cout << endl;
  }

  cout << "Time difference = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "[µs]" << endl;

  return 0;
}