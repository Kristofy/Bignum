#pragma once

#include <bits/stdc++.h>
#include "timer.h"

using namespace std;

class uint512_t_noavx {
    static const uint32_t LOG10_BASE = 9;
    static const uint32_t BASE = 1000000000;
    static constexpr uint32_t pow10[9] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};


public:
  uint512_t_noavx() = default;

  uint512_t_noavx(const int& n) {
    // TODO: check if n is valid
      digits.fill(0);
      digits[0] = n;
  }

  uint512_t_noavx(const string &s) {
    stringstream ss(s);
    ss >> *this;
  }

  uint512_t_noavx(const char *s) : uint512_t_noavx(string(s)) {}

  uint512_t_noavx(const uint512_t_noavx &other) = default;
  uint512_t_noavx(uint512_t_noavx &&other) = default;
  uint512_t_noavx &operator=(const uint512_t_noavx &other) = default;
  uint512_t_noavx &operator=(uint512_t_noavx &&other) = default;
  ~uint512_t_noavx() = default;

 
  inline auto operator+(const uint512_t_noavx &other) const -> uint512_t_noavx {
    plus++;
    plus_timer.resume();
    uint512_t_noavx result;
    uint32_t carry = 0;
    for (int i = 0; i < 16; i++) {
        result.digits[i] = digits[i] + other.digits[i] + carry;
        carry = result.digits[i] / BASE;
        result.digits[i] %= BASE;
    }
    plus_timer.stop();
    return result;
  }

  inline auto operator-(const uint512_t_noavx &other) const -> uint512_t_noavx {
    minus++;
    minus_timer.resume();
    uint512_t_noavx result;
    uint32_t carry = 0;
    for (int i = 0; i < 16; i++) {
        result.digits[i] = digits[i] - other.digits[i] - carry;
        carry = result.digits[i] / BASE;
        result.digits[i] %= BASE;
    }
    minus_timer.stop();
    return result;
  }

  inline auto operator*(const uint512_t_noavx &other) const -> uint512_t_noavx {
    mult++;
    mult_timer.resume();
    uint512_t_noavx result = 0;
    uint64_t carry[16 + 1];
    fill(carry, carry + 16 + 1, 0);
    for (int i = 0; i < 16; i++) {
      for (int j = 0; j <= i; j++) {
          uint64_t t = (uint64_t)digits[j] * other.digits[i - j];
          result.digits[i] += t % BASE;
          carry[i + 1] += t / BASE + result.digits[i] / BASE;
          result.digits[i] %= BASE;
      }
    }

    uint64_t c = 0;
    for (int i = 0; i < 16; i++) {
      uint64_t c2 = carry[i] + c;
      result.digits[i] += c2 % BASE;
      c = c2 / BASE + result.digits[i] / BASE;
      result.digits[i] %= BASE;
    }
    mult_timer.stop();
    return result;
  }

  inline auto operator/(const uint512_t_noavx &other) const -> uint512_t_noavx {
    div++;
    div_timer.resume();
    uint512_t_noavx result = 0;
    uint512_t_noavx remainder = 0;
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

  inline auto operator%(const uint512_t_noavx &other) const -> uint512_t_noavx {
    mod++;
    mod_timer.resume();
    uint512_t_noavx result = 0;
    uint512_t_noavx remainder = 0;
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

  inline auto operator==(const uint512_t_noavx &other) const -> bool {
    eq++;
    eq_timer.resume();
    for (int i = 0; i < 16; i++) {
      if (digits[i] != other.digits[i]) {
        return false;
      }
    }
    eq_timer.stop();
    return true;
  }

  inline auto operator!=(const uint512_t_noavx &other) const -> bool { return !(*this == other); }

  inline auto operator<(const uint512_t_noavx &other) const -> bool {
    lt_timer.resume();
    for (int i = 15; i >= 0; i--) {
      if (digits[i] != other.digits[i]) {
        return digits[i] < other.digits[i];
      }
    }
    lt_timer.stop();
    return false;
  }

  inline auto operator>(const uint512_t_noavx &other) const -> bool { return other < *this; }
  inline auto operator<=(const uint512_t_noavx &other) const -> bool { return !(other < *this); }
  inline auto operator>=(const uint512_t_noavx &other) const -> bool { return !(*this < other); }

private:
  array<uint32_t, 16> digits;
  
public:
  static uint64_t plus;
  static AdditiveTimer plus_timer;
  static uint64_t minus;
  static AdditiveTimer minus_timer;
  static uint64_t mult;
  static AdditiveTimer mult_timer;
  static uint64_t div;
  static AdditiveTimer div_timer;
  static uint64_t mod;
  static AdditiveTimer mod_timer;
  static uint64_t eq;
  static AdditiveTimer eq_timer;
  static uint64_t lt;
  static AdditiveTimer lt_timer;

  
  
  friend ostream &operator<<(ostream &os, const uint512_t_noavx &n) {
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

  friend istream &operator>>(istream &is, uint512_t_noavx &n) {
    string num_str; is >> num_str;
    
    n.digits.fill(0);

    // read in the 512 bits integer and store it in the array
    const int num_digits = num_str.size();
    // TODO: check if valid unsigned integer
    for(int i = 0; i < num_digits; i++){
        const int digit = num_str[i] - '0';
        const int digit_index = num_digits - i - 1;
        const int digit_array_index = digit_index / LOG10_BASE;
        const int digit_array_offset = digit_index % 9;
        n.digits[digit_array_index] += digit * pow10[digit_array_offset];
    }

    return is;
  }
};

uint64_t uint512_t_noavx::plus{};
AdditiveTimer uint512_t_noavx::plus_timer{};
uint64_t uint512_t_noavx::minus{};
AdditiveTimer uint512_t_noavx::minus_timer{};
uint64_t uint512_t_noavx::mult{};
AdditiveTimer uint512_t_noavx::mult_timer{};
uint64_t uint512_t_noavx::div{};
AdditiveTimer uint512_t_noavx::div_timer{};
uint64_t uint512_t_noavx::mod{};
AdditiveTimer uint512_t_noavx::mod_timer{};
uint64_t uint512_t_noavx::eq{};
AdditiveTimer uint512_t_noavx::eq_timer{};
uint64_t uint512_t_noavx::lt{};
AdditiveTimer uint512_t_noavx::lt_timer{};