#include <cstdint>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using namespace std;
using namespace std::chrono;

inline uint64_t my_popcnt(uint64_t x) {
  x = (x & 0x5555555555555555ULL) + ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1);
  x = (x & 0x3333333333333333ULL) + ((x & 0xCCCCCCCCCCCCCCCCULL) >> 2);
  x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x & 0xF0F0F0F0F0F0F0F0ULL) >> 4);
  x *= 0x0101010101010101ULL;
  return (x >> 56) & 0xFFULL;
}

int main() {
  random_device rnd_seed_gen;
  mt19937_64 rnd_engine(rnd_seed_gen());

  vector<uint64_t> rnd_values(1 << 16);
  for (uint64_t i = 0; i < rnd_values.size(); ++i) {
    rnd_values[i] = rnd_engine();
  }

  high_resolution_clock::time_point begin = high_resolution_clock::now();
  uint64_t n = 0;
  for (uint64_t m = 0; m < (1 << 30); m += rnd_values.size()) {
    for (uint64_t i = 0; i < rnd_values.size(); ++i) {
      n += __builtin_popcountll(rnd_values[i]);
    }
  }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  cout << "__builtin_popcountll: n = " << n
    << ", speed = " << ((1 << 30) / elapsed) << " values/ns" << endl;

  begin = high_resolution_clock::now();
  n = 0;
  for (uint64_t m = 0; m < (1 << 30); m += rnd_values.size()) {
    for (uint64_t i = 0; i < rnd_values.size(); ++i) {
      n += my_popcnt(rnd_values[i]);
    }
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  cout << "my_popcnt: n = " << n
    << ", speed = " << ((1 << 30) / elapsed) << " values/ns" << endl;
}
