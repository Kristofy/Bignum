#pragma once

#include <bits/stdc++.h>

using namespace std;

class uint512_t {
    static const uint32_t LOG10_BASE = 9;
    static const uint32_t BASE = 1000000000;
    static constexpr uint32_t pow10[9] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};


public:
  uint512_t() = default;

  uint512_t(const int& n) {
    // TODO: check if n is valid
      digits.fill(0);
      digits[0] = n;
  }

  uint512_t(const string &s) {
    stringstream ss(s);
    ss >> *this;
  }

  uint512_t(const char *s) : uint512_t(string(s)) {}

  uint512_t(const uint512_t &other) = default;
  uint512_t(uint512_t &&other) = default;
  uint512_t &operator=(const uint512_t &other) = default;
  uint512_t &operator=(uint512_t &&other) = default;
  ~uint512_t() = default;

  inline auto operator+(const uint512_t &other) const -> uint512_t {
    uint512_t result;
    uint32_t carry = 0;
    for (int i = 0; i < 16; i++) {
        result.digits[i] = digits[i] + other.digits[i] + carry;
        carry = result.digits[i] / BASE;
        result.digits[i] %= BASE;
    }
    return result;
  }

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
  array<uint32_t, 16> digits;

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