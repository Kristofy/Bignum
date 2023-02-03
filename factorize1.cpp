#include <bits/stdc++.h>
#include "timer.h"
using namespace std;

using bigint = long long;

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

auto factorize(bigint n) -> vector<bigint> {
  vector<bigint> factors;
  for (bigint i = 2; i * i <= n; i++) {
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
  vector<pair<bigint, vector<bigint>>> test_data = {{1223456789123LL, {}},   {1223456789126LL, {}},
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