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

  void stop() {
    stopped = true;
    auto end = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << name << " took " << elapsed << " ms" << std::endl;
  }

private:
    std::string name;
    std::chrono::steady_clock::time_point start;
    bool stopped = false;
};