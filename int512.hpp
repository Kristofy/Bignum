#pragma once

#pragma GCC target("avx2,bmi2")

#include <bits/stdc++.h>
#include <immintrin.h>

#include "timer.h"
using namespace std;


class FastDivider {
  private:
    enum {
        LIBDIVIDE_16_SHIFT_MASK = 0x1F,
        LIBDIVIDE_32_SHIFT_MASK = 0x1F,
        LIBDIVIDE_64_SHIFT_MASK = 0x3F,
        LIBDIVIDE_ADD_MARKER = 0x40,
        LIBDIVIDE_NEGATIVE_DIVISOR = 0x80
    };
  
    struct libdivide_u64_t {
        uint64_t magic;
        uint8_t more;
    };

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



    // libdivide_128_div_64_to_64: divides a 128-bit uint {numhi, numlo} by a 64-bit uint {den}. The
    // result must fit in 64 bits. Returns the quotient directly and the remainder in *r
    static inline uint64_t libdivide_128_div_64_to_64(
        uint64_t numhi, uint64_t numlo, uint64_t den, uint64_t *r) {
        // N.B. resist the temptation to use __uint128_t here.
        // In LLVM compiler-rt, it performs a 128/128 -> 128 division which is many times slower than
        // necessary. In gcc it's better but still slower than the divlu implementation, perhaps because
        // it's not LIBDIVIDE_INLINEd.
    #if defined(LIBDIVIDE_X86_64) && defined(LIBDIVIDE_GCC_STYLE_ASM)
        uint64_t result;
        __asm__("divq %[v]" : "=a"(result), "=d"(*r) : [v] "r"(den), "a"(numlo), "d"(numhi));
        return result;
    #else
        // We work in base 2**32.
        // A uint32 holds a single digit. A uint64 holds two digits.
        // Our numerator is conceptually [num3, num2, num1, num0].
        // Our denominator is [den1, den0].
        const uint64_t b = ((uint64_t)1 << 32);

        // The high and low digits of our computed quotient.
        uint32_t q1;
        uint32_t q0;

        // The normalization shift factor.
        int shift;

        // The high and low digits of our denominator (after normalizing).
        // Also the low 2 digits of our numerator (after normalizing).
        uint32_t den1;
        uint32_t den0;
        uint32_t num1;
        uint32_t num0;

        // A partial remainder.
        uint64_t rem;

        // The estimated quotient, and its corresponding remainder (unrelated to true remainder).
        uint64_t qhat;
        uint64_t rhat;

        // Variables used to correct the estimated quotient.
        uint64_t c1;
        uint64_t c2;

        // Check for overflow and divide by 0.
        if (numhi >= den) {
            if (r != NULL) *r = ~0ull;
            return ~0ull;
        }

        // Determine the normalization factor. We multiply den by this, so that its leading digit is at
        // least half b. In binary this means just shifting left by the number of leading zeros, so that
        // there's a 1 in the MSB.
        // We also shift numer by the same amount. This cannot overflow because numhi < den.
        // The expression (-shift & 63) is the same as (64 - shift), except it avoids the UB of shifting
        // by 64. The funny bitwise 'and' ensures that numlo does not get shifted into numhi if shift is
        // 0. clang 11 has an x86 codegen bug here: see LLVM bug 50118. The sequence below avoids it.
        shift = __builtin_clzll(den);
        den <<= shift;
        numhi <<= shift;
        numhi |= (numlo >> (-shift & 63)) & (-(int64_t)shift >> 63);
        numlo <<= shift;

        // Extract the low digits of the numerator and both digits of the denominator.
        num1 = (uint32_t)(numlo >> 32);
        num0 = (uint32_t)(numlo & 0xFFFFFFFFu);
        den1 = (uint32_t)(den >> 32);
        den0 = (uint32_t)(den & 0xFFFFFFFFu);

        // We wish to compute q1 = [n3 n2 n1] / [d1 d0].
        // Estimate q1 as [n3 n2] / [d1], and then correct it.
        // Note while qhat may be 2 digits, q1 is always 1 digit.
        qhat = numhi / den1;
        rhat = numhi % den1;
        c1 = qhat * den0;
        c2 = rhat * b + num1;
        if (c1 > c2) qhat -= (c1 - c2 > den) ? 2 : 1;
        q1 = (uint32_t)qhat;

        // Compute the true (partial) remainder.
        rem = numhi * b + num1 - q1 * den;

        // We wish to compute q0 = [rem1 rem0 n0] / [d1 d0].
        // Estimate q0 as [rem1 rem0] / [d1] and correct it.
        qhat = rem / den1;
        rhat = rem % den1;
        c1 = qhat * den0;
        c2 = rhat * b + num0;
        if (c1 > c2) qhat -= (c1 - c2 > den) ? 2 : 1;
        q0 = (uint32_t)qhat;

        // Return remainder if requested.
        if (r != NULL) *r = (rem * b + num0 - q0 * den) >> shift;
        return ((uint64_t)q1 << 32) | q0;
    #endif
    }

    static inline struct libdivide_u64_t libdivide_internal_u64_gen(uint64_t d, int branchfree) {
      

        struct libdivide_u64_t result;
        uint32_t floor_log_2_d = 63 - __builtin_clzll(d);

        // Power of 2
        if ((d & (d - 1)) == 0) {
            // We need to subtract 1 from the shift value in case of an unsigned
            // branchfree divider because there is a hardcoded right shift by 1
            // in its division algorithm. Because of this we also need to add back
            // 1 in its recovery algorithm.
            result.magic = 0;
            result.more = (uint8_t)(floor_log_2_d - (branchfree != 0));
        } else {
            uint64_t proposed_m, rem;
            uint8_t more;
            // (1 << (64 + floor_log_2_d)) / d
            proposed_m = libdivide_128_div_64_to_64((uint64_t)1 << floor_log_2_d, 0, d, &rem);

            const uint64_t e = d - rem;

            // This power works if e < 2**floor_log_2_d.
            if (!branchfree && e < ((uint64_t)1 << floor_log_2_d)) {
                // This power works
                more = (uint8_t)floor_log_2_d;
            } else {
                // We have to use the general 65-bit algorithm.  We need to compute
                // (2**power) / d. However, we already have (2**(power-1))/d and
                // its remainder. By doubling both, and then correcting the
                // remainder, we can compute the larger division.
                // don't care about overflow here - in fact, we expect it
                proposed_m += proposed_m;
                const uint64_t twice_rem = rem + rem;
                if (twice_rem >= d || twice_rem < rem) proposed_m += 1;
                more = (uint8_t)(floor_log_2_d | LIBDIVIDE_ADD_MARKER);
            }
            result.magic = 1 + proposed_m;
            result.more = more;
            // result.more's shift should in general be ceil_log_2_d. But if we
            // used the smaller power, we subtract one from the shift because we're
            // using the smaller power. If we're using the larger power, we
            // subtract one from the shift because it's taken care of by the add
            // indicator. So floor_log_2_d happens to be correct in both cases,
            // which is why we do it outside of the if statement.
        }
        return result;
    }


    struct libdivide_u64_t libdivide_u64_gen(uint64_t d) {
        return libdivide_internal_u64_gen(d, 0);
    }

    static inline __m256i libdivide_u64_do_vec256(__m256i numers, const struct libdivide_u64_t *denom) {
        uint8_t more = denom->more;
        if (!denom->magic) {
            return _mm256_srli_epi64(numers, more);
        } else {
            __m256i q = libdivide_mullhi_u64_vec256(numers, _mm256_set1_epi64x(denom->magic));
            if (more & LIBDIVIDE_ADD_MARKER) {
                // uint32_t t = ((numer - q) >> 1) + q;
                // return t >> denom->shift;
                uint32_t shift = more & LIBDIVIDE_64_SHIFT_MASK;
                __m256i t = _mm256_add_epi64(_mm256_srli_epi64(_mm256_sub_epi64(numers, q), 1), q);
                return _mm256_srli_epi64(t, shift);
            } else {
                return _mm256_srli_epi64(q, more);
            }
        }
    }


    public:
      FastDivider(uint64_t denom)
        : denom(libdivide_u64_gen(denom)) 
      {}

      inline __m256i div_x8(__m256i num) const {
        return libdivide_u64_do_vec256(num, &denom);
      }

    private:
      libdivide_u64_t denom;

    
};



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
    plus_timer.resume();
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

    plus_timer.stop();
    return result;
  }

  // TODO: Optimize
  inline auto operator-(const uint512_t &other) const -> uint512_t {
    uint512_t result;
    minus_timer.resume();
    uint32_t carry = 0;
    for (int i = 0; i < 16; i++) {
        result.digits[i] = digits[i] - other.digits[i] - carry;
        carry = result.digits[i] / BASE;
        result.digits[i] %= BASE;
    }

    minus_timer.stop();
    return result;
  }

  inline auto operator*(const uint512_t &other) const -> uint512_t {
    uint512_t result = 0;
    mult_timer.resume();
    alignas(32) static uint64_t carry[16];

    static const __m256i ZeroX8 = _mm256_set1_epi32(0);
    static const __m256i Base01X4 = _mm256_set_epi32(0, BASE, 0, BASE, 0, BASE, 0, BASE);
    static const __m256i MaxIndexX8 = _mm256_set1_epi32(15);
    static const __m256i ShiftRMask = _mm256_set_epi32(6,5,4,3,2,1,0,7);
    static const FastDivider BaseDivider(BASE);
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
      C = _mm256_add_epi64(A, B);

      carry[i] = (uint64_t)_mm256_extract_epi64(C, 0) + (uint64_t)_mm256_extract_epi64(C, 2);
    }

    A = _mm256_load_si256((__m256i*)(carry + 0));
    B = _mm256_load_si256((__m256i*)(carry + 4));
    C = _mm256_load_si256((__m256i*)(carry + 8));
    D = _mm256_load_si256((__m256i*)(carry + 12));

    X = BaseDivider.div_x8(A);
    Y = BaseDivider.div_x8(B);
    Z = BaseDivider.div_x8(C);
    W = BaseDivider.div_x8(D);

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
    uint64_t q;
    for (int i = 1; i < 16; i++) {
      q = c / BASE;
      result.digits[i] += c - q * BASE;
      c = q + (result.digits[i] >= BASE);
      if(result.digits[i] >= BASE) result.digits[i] -= BASE;
      c += carry[i];
    }


    mult_timer.stop();
    return result;
  }

  inline auto operator/(const uint512_t &other) const -> uint512_t {
    uint512_t result = 0;
    div_timer.resume();
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

    div_timer.stop();
    return result;
  }

  inline auto operator%(const uint512_t &other) const -> uint512_t {
    uint512_t result = 0;
    mod_timer.resume();
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

    mod_timer.stop();
    return remainder;
  }

  inline auto operator==(const uint512_t &other) const -> bool {
    eq_timer.resume();
    for (int i = 0; i < 16; i++) {
      if (digits[i] != other.digits[i]) {
        return false;
      }
    }
    eq_timer.stop();
    return true;
  }

  inline auto operator!=(const uint512_t &other) const -> bool { return !(*this == other); }

  inline auto operator<(const uint512_t &other) const -> bool {
    lt_timer.resume();
    for (int i = 15; i >= 0; i--) {
      if (digits[i] != other.digits[i]) {
        return digits[i] < other.digits[i];
      }
    }
    lt_timer.stop();
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

  static AdditiveTimer plus_timer;
  static AdditiveTimer minus_timer;
  static AdditiveTimer mult_timer;
  static AdditiveTimer div_timer;
  static AdditiveTimer mod_timer;
  static AdditiveTimer eq_timer;
  static AdditiveTimer lt_timer;
};

AdditiveTimer uint512_t::plus_timer{};
AdditiveTimer uint512_t::minus_timer{};
AdditiveTimer uint512_t::mult_timer{};
AdditiveTimer uint512_t::div_timer{};
AdditiveTimer uint512_t::mod_timer{};
AdditiveTimer uint512_t::eq_timer{};
AdditiveTimer uint512_t::lt_timer{};