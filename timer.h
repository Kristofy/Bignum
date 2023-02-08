#pragma once

#include <chrono>
#include <iostream>
#include <string>
#include <vector>



class Timer {
public:
  Timer(const std::string &name) : name(name) {
    start = std::chrono::steady_clock::now();
  }

  ~Timer() {
    if (!stopped) {
      stop();
    }
  }

  auto stop() -> int64_t {
    stopped = true;
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << name << " took " << elapsed/1000 << " ms" << std::endl;
    return elapsed;
  }

private:
    std::string name;
    std::chrono::steady_clock::time_point start;
    bool stopped = false;
};

class AdditiveTimer {
public:
  void resume(){
    start = std::chrono::steady_clock::now();
  }

  auto stop() {
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    total += elapsed;
    return elapsed;
  }

  auto get_total() -> int64_t {
    return total;
  }

  private:
    std::chrono::steady_clock::time_point start;
    uint64_t total = 0;
};
