#include <bits/stdc++.h>

#include "int512.hpp"
#include "int512_noavx.hpp"
#include "timer.h"

using namespace std;

#define DONT_OPTIMIZE(x)                                                                                               \
  { void *volatile dummy = &x; }

template <typename T> auto fib_test() -> int64_t {
  Timer t("test_addition_fib");
  for (int j = 0; j < 100000; j++) {
    T a = 1, b = 1, c;
    DONT_OPTIMIZE(a);
    DONT_OPTIMIZE(b);
    DONT_OPTIMIZE(c);
    // fibonacci numbers
    for (int i = 0; i < 650; i++) {
      c = a + b;
      a = b;
      b = c;
    }
  }
  return t.stop();
}

template <typename T> auto gauss_test() -> int64_t {
  Timer t("test_addition_gauss");
  T sum = 0;
  DONT_OPTIMIZE(sum);
  const int n = 100000000;
  for (int i = 1; i <= n; i++) {
    sum = sum + i;
  }
  if (sum + sum != T(n) * T(n + 1)) {
    cout << "wrong result" << endl;
  }

  return t.stop();
}

template <typename T> auto factorial_test() -> int64_t {
  Timer t("test_multiplication_factorial");
  for (int j = 0; j < 50000; j++) {
    T prod = 1;
    DONT_OPTIMIZE(prod);
    const int n = 90;
    for (int i = 1; i <= n; i++) {
      prod = prod * i;
    }
  }

  return t.stop();
}

template <typename T> auto collatz_test() -> int64_t {
  Timer t("test_division_collatz");
  T n = "338838125384605298683130670138635876715138828490733673519118221753294017425139038491279";
  DONT_OPTIMIZE(n);
  T two = 2;
  while (n != 1) {
    if (n % two == 0) {
      n = n / two;
    } else {
      n = n * 3 + 1;
    }
  }

  return t.stop();
}

template <typename T> auto mersenne_prime_test() -> int64_t {
  Timer t("test_division_mersenne_prime");
  // Find the first n mersenne prime
  const int n = 8;
  auto is_prime = [](const T &p) -> bool {
    if (p == 1)
      return false;
    if (p == 2)
      return true;
    if (p % 2 == 0)
      return false;
    T i = 3;
    while (i * i <= p) {
      if (p - p / i * i == 0)
        return false;
      i = i + 2;
    }
    return true;
  };

  T p = 3;
  T two_to_i = 2;
  int twop = 1;
  DONT_OPTIMIZE(p);
  for (int i = 1; i <= n;) {
    two_to_i = two_to_i * 2;
    twop++;
    p = two_to_i - 1;
    if (is_prime(p)) {
      i++;
    }
  }
  return t.stop();
}

int main() {

  cout << "Without avx:" << endl;
  auto noavx_fib = fib_test<uint512_t_noavx>();
  auto noavx_gauss = gauss_test<uint512_t_noavx>();
  auto noavx_factorial = factorial_test<uint512_t_noavx>();
  auto noavx_collatz = collatz_test<uint512_t_noavx>();
  auto noavx_mersenne_prime = mersenne_prime_test<uint512_t_noavx>();
  cout << endl;

  cout << "With avx:" << endl;
  auto avx_fib = fib_test<uint512_t>();
  auto avx_gauss = gauss_test<uint512_t>();
  auto avx_factorial = factorial_test<uint512_t>();
  auto avx_collatz = collatz_test<uint512_t>();
  auto avx_mersenne_prime = mersenne_prime_test<uint512_t>();
  cout << endl;

  cout << "Comparisons:" << endl;
  cout << "fib: \t" << avx_fib/1000.0 << "\t" << noavx_fib/1000.0 << "\t" << (double)avx_fib / noavx_fib << endl;
  cout << "gauss: \t" << avx_gauss/1000.0 << "\t" << noavx_gauss/1000.0 << "\t" << (double)avx_gauss / noavx_gauss << endl;
  cout << "factorial: \t" << avx_factorial/1000.0 << "\t" << noavx_factorial/1000.0 << "\t" << (double)avx_factorial / noavx_factorial << endl;
  cout << "collatz: \t" << avx_collatz/1000.0 << "\t" << noavx_collatz/1000.0 << "\t" << (double)avx_collatz / noavx_collatz << endl;
  cout << "mersenne_prime: \t" << avx_mersenne_prime/1000.0 << "\t" << noavx_mersenne_prime/1000.0 << "\t" << (double)avx_mersenne_prime / noavx_mersenne_prime << endl;
  cout << endl;


  cout << "Operators:" << endl;
  cout << "plus: \t" << uint512_t::plus_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::plus_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::plus_timer.get_total() / (double)uint512_t::plus_timer.get_total() << endl;
  cout << "minus: \t" << uint512_t::minus_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::minus_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::minus_timer.get_total() / (double)uint512_t::minus_timer.get_total() << endl;
  cout << "mult: \t" << uint512_t::mult_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::mult_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::mult_timer.get_total() / (double)uint512_t::mult_timer.get_total() << endl;
  cout << "div: \t" << uint512_t::div_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::div_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::div_timer.get_total() / (double)uint512_t::div_timer.get_total() << endl;
  cout << "mod: \t" << uint512_t::mod_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::mod_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::mod_timer.get_total() / (double)uint512_t::mod_timer.get_total() << endl;
  cout << "eq: \t" << uint512_t::eq_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::eq_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::eq_timer.get_total() / (double)uint512_t::eq_timer.get_total() << endl;
  cout << "lt: \t" << uint512_t::lt_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::lt_timer.get_total() / 1000.0 << "\t" << uint512_t_noavx::lt_timer.get_total() / (double)uint512_t::lt_timer.get_total() << endl;
  cout << endl;



  return 0;
}