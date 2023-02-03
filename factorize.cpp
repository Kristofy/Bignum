#pragma GCC target("avx2,bmi2")

#include <bits/stdc++.h>
#include <immintrin.h>

#include "timer.h"
using namespace std;


// Here, y is assumed to contain one 64-bit value repeated.
static inline __m256i libdivide_mullhi_u64_vec256(__m256i x, __m256i y) {
    // see m128i variant for comments.
    __m256i x0y0 = _mm256_mul_epu32(x, y);
    __m256i x0y0_hi = _mm256_srli_epi64(x0y0, 32);

    __m256i x1 = _mm256_shuffle_epi32(x, _MM_SHUFFLE(3, 3, 1, 1));
    __m256i y1 = _mm256_shuffle_epi32(y, _MM_SHUFFLE(3, 3, 1, 1));

    __m256i x0y1 = _mm256_mul_epu32(x, y1);
    __m256i x1y0 = _mm256_mul_epu32(x1, y);
    __m256i x1y1 = _mm256_mul_epu32(x1, y1);

    __m256i mask = _mm256_set1_epi64x(0xFFFFFFFF);
    __m256i temp = _mm256_add_epi64(x1y0, x0y0_hi);
    __m256i temp_lo = _mm256_and_si256(temp, mask);
    __m256i temp_hi = _mm256_srli_epi64(temp, 32);

    temp_lo = _mm256_srli_epi64(_mm256_add_epi64(temp_lo, x0y1), 32);
    temp_hi = _mm256_add_epi64(x1y1, temp_hi);
    return _mm256_add_epi64(temp_lo, temp_hi);
}

static inline uint64_t libdivide_128_div_64_to_64(uint64_t numhi, uint64_t numlo, uint64_t den, uint64_t *r) {
    uint64_t result;
    __asm__("divq %[v]" : "=a"(result), "=d"(*r) : [v] "r"(den), "a"(numlo), "d"(numhi));
    return result;
}

static inline uint64_t libdivide_u64_gen(uint64_t d) {
    uint64_t result;
    uint64_t rem;
    result = 1 + 2*libdivide_128_div_64_to_64((uint64_t)1 << 29, 0, d, &rem);
    const uint64_t twice_rem = rem + rem;
    if (twice_rem >= d || twice_rem < rem) result += 1;
    return result;
}


__m256i libdivide_u64_do_vec256(__m256i numers, uint64_t denom) {
    __m256i q = libdivide_mullhi_u64_vec256(numers, _mm256_set1_epi64x(denom));
    __m256i t = _mm256_add_epi64(_mm256_srli_epi64(_mm256_sub_epi64(numers, q), 1), q);
    return _mm256_srli_epi64(t, 29);
}

// 522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139

// it stores the number in an array of 16 uint32_t each up to 9 digits
// the total maximum is 16 * 9 = 144 digits
class uint512_t {
    static const uint32_t LOG10_BASE = 9;
    static const uint32_t BASE = 1000000000;
    static constexpr uint32_t pow10[9] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};


public:
  uint512_t() = default;

  uint512_t(const int& n) {
    // TODO: check if n is valid
      fill(digits, digits + 16, 0);
      digits[0] = n;
  }

  uint512_t(const string &s) {
    fill(digits, digits + 16, 0);

    // read in the 512 bits integer and store it in the array
    const int num_digits = s.size();
    // TODO: check if valid unsigned integer
    for(int i = 0; i < num_digits; i++){
        const int digit = s[i] - '0';
        const int digit_index = num_digits - i - 1;
        const int digit_array_index = digit_index / LOG10_BASE;
        const int digit_array_offset = digit_index % 9;
        digits[digit_array_index] += digit * pow10[digit_array_offset];
    }

  }

  uint512_t(const char *s) : uint512_t(string(s)) {}

  uint512_t(const uint512_t &other) = default;
  uint512_t(uint512_t &&other) = default;
  uint512_t &operator=(const uint512_t &other) = default;
  uint512_t &operator=(uint512_t &&other) = default;
  ~uint512_t() = default;

  auto operator+(const uint512_t &other) const -> uint512_t {
    uint512_t result;
    static const __m256i BaseX8 = _mm256_set1_epi32(BASE);
    static const __m256i OneX8 = _mm256_set1_epi32(1);
    static const __m256i ShiftRMask = _mm256_set_epi32(6,5,4,3,2,1,0,7);

    __m256i A, B, C, D;
    __m256i X, Y;
    __m256i CL, CH;

    X = _mm256_add_epi32(Lo, other.Lo);
    Y = _mm256_add_epi32(Hi, other.Hi);

    CL = _mm256_cmpgt_epi32(BaseX8, X);
    CH = _mm256_cmpgt_epi32(BaseX8, Y);

    A = _mm256_sub_epi32(X, BaseX8);
    B = _mm256_and_si256(CL, X);
    C = _mm256_andnot_si256(CL, A);
    X = _mm256_or_si256(B, C);

    A = _mm256_sub_epi32(Y, BaseX8);
    B = _mm256_and_si256(CH, Y);
    C = _mm256_andnot_si256(CH, A);
    Y = _mm256_or_si256(B, C);
  
    A = _mm256_andnot_si256(CL, OneX8);
    B = _mm256_andnot_si256(CH, OneX8);

    B = _mm256_insert_epi32(B, _mm256_extract_epi32(A, 7), 7);
    A = _mm256_insert_epi32(A, 0, 7);
    CL = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(A), ShiftRMask));
    CH = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(B), ShiftRMask));

    result.Lo = _mm256_add_epi32(X, CL);
    result.Hi = _mm256_add_epi32(Y, CH);
    
    CL = _mm256_cmpgt_epi32(BaseX8, X);
    CH = _mm256_cmpgt_epi32(BaseX8, Y);

    A = _mm256_or_si256(CL, CH);
    bool is_carry =  _mm256_testz_si256(A, A) == 0;

    if(is_carry){
      uint32_t c = 0;
      for (int i = 0; c && i < 16; i++) {
          result.digits[i] += c;
          c = result.digits[i] >= BASE;
          if(c) result.digits[i] -= BASE;
      }
    }

    return result;
  }

  // TODO: Optimize
  inline auto operator-(const uint512_t &other) const -> uint512_t {
    uint512_t result;
    uint32_t carry = 0;
    for (int i = 0; i < 16; i++) {
        result.digits[i] = digits[i] - other.digits[i] - carry;
        carry = result.digits[i] / BASE;
        result.digits[i] %= BASE;
    }
    return result;
  }

  inline auto operator*(const uint512_t &other) const -> uint512_t {
    uint512_t result = 0;
    alignas(64) static uint64_t carry[16];

    static const __m256i ZeroX8 = _mm256_set1_epi32(0);
    static const __m256i Base01X4 = _mm256_set_epi32(0, BASE, 0, BASE, 0, BASE, 0, BASE);
    static const __m256i MaxIndexX8 = _mm256_set1_epi32(15);
    static const __m256i ShiftRMask = _mm256_set_epi32(6,5,4,3,2,1,0,7);
    static const uint64_t BaseDivider = libdivide_u64_gen(BASE);
    static const __m256i MaskGeneratorLo = _mm256_set_epi32(6, 5, 4, 3, 2, 1, 0, -1);
    static const __m256i MaskGeneratorHi = _mm256_set_epi32(14, 13, 12, 11, 10, 9, 8, 7);
    static const __m256i Mask01X8 = _mm256_set_epi32(0, -1, 0, -1, 0, -1, 0, -1);
    static const __m256i Mask01X4 = _mm256_set_epi32(0, 0, -1, -1, 0, 0, -1, -1);
    static const __m256i MultiplyShuffleMask = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
    static const __m256i ReverseShuffle = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    static const __m256i UnpackLowShuffle = _mm256_set_epi32(7, 5, 3, 1, 6, 4, 2, 0);
    static const __m256i UnpackHiShuffle = _mm256_set_epi32(6, 4, 2, 0, 7, 5, 3, 1);
    static const __m256i JOffsetLo = _mm256_set_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    static const __m256i JOffsetHi = _mm256_set_epi32(8, 9, 10, 11, 12, 13, 14, 15);

    union {
      alignas(64) uint64_t d[4];
      __m256i X;
    } R;

    __m256i A, B, C, D, BR, DR;
    __m256i X, Y, Z, W;
    __m256i ML0, ML1, MH0, MH1;
    __m256i MaskLo, MaskHi;
    __m256i IX8;
    BR = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(other.Hi), ReverseShuffle));
    DR = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(other.Lo), ReverseShuffle));

    for (int i = 0; i < 16; i++) {
      IX8 = _mm256_set1_epi32(i);
      MaskLo = _mm256_cmpgt_epi32(IX8, MaskGeneratorLo);
      MaskHi = _mm256_cmpgt_epi32(IX8, MaskGeneratorHi);


      DR = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(DR), ShiftRMask));
      int c = _mm256_extract_epi32(DR, 0);
      BR = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(BR), ShiftRMask));
      DR = _mm256_insert_epi32(DR, _mm256_extract_epi32(BR, 0) , 0);
      BR = _mm256_insert_epi32(BR, c, 0);


      A = _mm256_and_si256(this->Lo, MaskLo);
      B = _mm256_and_si256(BR, MaskLo);
      C = _mm256_and_si256(this->Hi, MaskHi);
      D = _mm256_and_si256(DR, MaskHi);
     
      // Set up for multiplication
      A = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(A), MultiplyShuffleMask));
      B = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(B), MultiplyShuffleMask));
      C = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(C), MultiplyShuffleMask));
      D = _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(D), MultiplyShuffleMask));
        
      // multiply the lower halfes
      X = _mm256_and_si256(A, Mask01X8);
      Y = _mm256_srli_epi64(_mm256_andnot_si256(Mask01X8, A), 32);
      Z = _mm256_and_si256(B, Mask01X8);
      W = _mm256_srli_epi64(_mm256_andnot_si256(Mask01X8, B), 32);

      ML0 = _mm256_mul_epu32(X, Z);
      ML1 = _mm256_mul_epu32(Y, W);

      // multiply the upper halfes
      X = _mm256_and_si256(C, Mask01X8);
      Y = _mm256_srli_epi64(_mm256_andnot_si256(Mask01X8, C), 32);
      Z = _mm256_and_si256(D, Mask01X8);
      W = _mm256_srli_epi64(_mm256_andnot_si256(Mask01X8, D), 32);

      MH0 = _mm256_mul_epu32(X, Z);
      MH1 = _mm256_mul_epu32(Y, W);

      A = _mm256_add_epi64(ML0, ML1);
      B = _mm256_add_epi64(MH0, MH1);
      C = _mm256_add_epi64(A, B);

      A = _mm256_and_si256(Mask01X4, C);
      B = _mm256_srli_si256(_mm256_andnot_si256(Mask01X4, C), 8);
      R.X = _mm256_add_epi64(A, B);

      // A = _mm256_permute2f128_si256(C, C, 0x00);
      // B = _mm256_permute2f128_si256(C, C, 0x11);
      // carry[i] = _mm256_extract_epi64(_mm256_add_epi64(A, B), 0);
      carry[i] = R.d[0] + R.d[2];
    }

    A = _mm256_load_si256((__m256i*)(carry + 0));
    B = _mm256_load_si256((__m256i*)(carry + 4));
    C = _mm256_load_si256((__m256i*)(carry + 8));
    D = _mm256_load_si256((__m256i*)(carry + 12));

    X = libdivide_u64_do_vec256(A, BaseDivider);
    Y = libdivide_u64_do_vec256(B, BaseDivider);
    Z = libdivide_u64_do_vec256(C, BaseDivider);
    W = libdivide_u64_do_vec256(D, BaseDivider);

    ML0 = _mm256_mul_epu32(X, Base01X4);
    ML1 = _mm256_mul_epu32(Y, Base01X4);
    MH0 = _mm256_mul_epu32(Z, Base01X4);
    MH1 = _mm256_mul_epu32(W, Base01X4);

    A = _mm256_sub_epi64(A, ML0);
    B = _mm256_sub_epi64(B, ML1);
    C = _mm256_sub_epi64(C, MH0);
    D = _mm256_sub_epi64(D, MH1);

    result.Lo = _mm256_or_si256(
        _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(A), UnpackLowShuffle)),
        _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(B), UnpackHiShuffle)));
    result.Hi = _mm256_or_si256(
        _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(C), UnpackLowShuffle)),
        _mm256_castps_si256(_mm256_permutevar8x32_ps(_mm256_castsi256_ps(D), UnpackHiShuffle)));


    _mm256_store_si256((__m256i*)(carry + 0), X);
    _mm256_store_si256((__m256i*)(carry + 4), Y);
    _mm256_store_si256((__m256i*)(carry + 8), Z);
    _mm256_store_si256((__m256i*)(carry + 12), W);

    uint64_t c = carry[0];
    for (int i = 1; i < 16; i++) {
      result.digits[i] += c % BASE;
      c = (c / BASE) + (result.digits[i] >= BASE);
      if(result.digits[i] >= BASE) result.digits[i] -= BASE;
      c += carry[i];
    }

    return result;
  }

  inline auto operator/(const uint512_t &other) const -> uint512_t {
    uint512_t result = 0;
    uint512_t remainder = 0;
    for (int i = 15; i >= 0; i--) {
      remainder = remainder * BASE + digits[i];
      uint32_t low = 0, high = BASE;
      while (low < high) {
        uint32_t mid = (low + high + 1) / 2;
        if (other * mid <= remainder) {
          low = mid;
        } else {
          high = mid - 1;
        }
      }
      result.digits[i] = low;
      remainder = remainder - other * low;
    }
    return result;
  }

  inline auto operator%(const uint512_t &other) const -> uint512_t {
    uint512_t result = 0;
    uint512_t remainder = 0;
    for (int i = 15; i >= 0; i--) {
      remainder = remainder * BASE + digits[i];
      uint32_t low = 0, high = BASE;
      while (low < high) {
        uint32_t mid = (low + high + 1) / 2;
        if (other * mid <= remainder) {
          low = mid;
        } else {
          high = mid - 1;
        }
      }
      result.digits[i] = low;
      remainder = remainder - other * low;
    }
    return remainder;
  }

  inline auto operator==(const uint512_t &other) const -> bool {
    for (int i = 0; i < 16; i++) {
      if (digits[i] != other.digits[i]) {
        return false;
      }
    }
    return true;
  }

  inline auto operator!=(const uint512_t &other) const -> bool { return !(*this == other); }

  inline auto operator<(const uint512_t &other) const -> bool {
    for (int i = 15; i >= 0; i--) {
      if (digits[i] != other.digits[i]) {
        return digits[i] < other.digits[i];
      }
    }
    return false;
  }

  inline auto operator>(const uint512_t &other) const -> bool { return other < *this; }
  inline auto operator<=(const uint512_t &other) const -> bool { return !(other < *this); }
  inline auto operator>=(const uint512_t &other) const -> bool { return !(*this < other); }

private:
  union {
     alignas(32) uint32_t digits[16];
     struct {
       __m256i Lo;
       __m256i Hi;
     };
   
  };
public:
  
  
  friend ostream &operator<<(ostream &os, const uint512_t &n) {
    // write to os without leading zeros
    bool first = true;
    for (int i = 15; i >= 0; i--) {
      if (first) {
        if (n.digits[i] != 0) {
          os << n.digits[i];
          first = false;
        }
      } else {
        os << setw(9) << setfill('0') << n.digits[i];
      }
    }
    return os;
  }

  friend istream &operator>>(istream &is, uint512_t &n) {
    string num_str; is >> num_str;
    n = uint512_t(num_str);
    return is;
  }
};

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
  uint512_t p = "522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139";
  int c = 0;
  for (uint512_t i = 3; i * i <= p; i = i + 2) {
    while (p % i == 0) {
      p = p / i;
      cout << "Found factor " << i << endl;
    }

    c++;
    if (c == 10000) {
      cout << "Checked up until " << i << " the ramaining number is: " << p << endl;
      c = 0;
    }
  }

  cout << "The biggest factor of p is " << p << endl;

  // fib_test();
  // gauss_test();
  // factorial_test();
  // collatz_test();
  // mersenne_prime_test();
  return 0;
}