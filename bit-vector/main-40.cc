#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <x86intrin.h>

using namespace std;
using namespace std::chrono;

class BitVector {
 public:
  BitVector()
    : words_(), n_bits_(0),
      ranks_(), select0s_(), select1s_(), n_zeros_(0), n_ones_(0) {}
  ~BitVector() {}

  uint64_t get(uint64_t i) const {
    assert(i < n_bits_);
    return (words_[i / 64] >> (i % 64)) & 1;
  }
  void set(uint64_t i, uint64_t bit) {
    assert(i < n_bits_);
    set_(i, bit);
  }

  void add(uint64_t bit) {
    if (n_bits_ % 256 == 0) {
      words_.resize((n_bits_ + 256) / 64);
    }
    set_(n_bits_, bit);
    ++n_bits_;
  }

  uint64_t n_bits() const {
    return n_bits_;
  }

  void build() {
    uint64_t n_blocks = words_.size() / 4;
    ranks_.resize(n_blocks + 1);
    uint64_t n_zeros = 0;
    uint64_t n_ones = 0;
    for (uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      ranks_[block_id].set_abs(n_ones);
      for (uint64_t j = 0; j < 4; ++j) {
        if (j != 0) {
          uint64_t rel = n_ones - ranks_[block_id].abs();
          ranks_[block_id].rels[j - 1] = rel;
        }

        uint64_t word_id = (block_id * 4) + j;
        uint64_t n_pops = __builtin_popcountll(words_[word_id]);
        uint64_t new_n_zeros = n_zeros + 64 - n_pops;
        if (((n_zeros + 255) / 256) != ((new_n_zeros + 255) / 256)) {
          uint64_t count = n_zeros;
          uint64_t word = ~words_[word_id];
          while (word != 0) {
            uint64_t pos = __builtin_ctzll(word);
            if (count % 256 == 0) {
              select0s_.push_back(((word_id * 64) + pos) >> 8);
              break;
            }
            word ^= 1UL << pos;
            ++count;
          }
        }
        n_zeros = new_n_zeros;

        uint64_t new_n_ones = n_ones + n_pops;
        if (((n_ones + 255) / 256) != ((new_n_ones + 255) / 256)) {
          uint64_t count = n_ones;
          uint64_t word = words_[word_id];
          while (word != 0) {
            uint64_t pos = __builtin_ctzll(word);
            if (count % 256 == 0) {
              select1s_.push_back(((word_id * 64) + pos) >> 8);
              break;
            }
            word ^= 1UL << pos;
            ++count;
          }
        }
        n_ones = new_n_ones;
      }
    }
    ranks_.back().set_abs(n_ones);
    select0s_.push_back(words_.size() * 64 / 256);
    select1s_.push_back(words_.size() * 64 / 256);
  }

  uint64_t rank(uint64_t i) const {
    uint64_t word_id = i / 64;
    uint64_t bit_id = i % 64;
    uint64_t rank_id = word_id / 4;
    uint64_t rel_id = word_id % 4;
    uint64_t n = ranks_[rank_id].abs();
    if (rel_id != 0) {
      n += ranks_[rank_id].rels[rel_id - 1];
    }
    n += __builtin_popcountll(words_[word_id] & ((1UL << bit_id) - 1));
    return n;
  }

  uint64_t select0(uint64_t i) const {
    const uint64_t block_id = i / 256;
    uint64_t begin = select0s_[block_id];
    uint64_t end = select0s_[block_id + 1] + 1;
    if (begin + 10 >= end) {
      while (i >= (256 * (begin + 1)) - ranks_[begin + 1].abs()) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < (256 * middle) - ranks_[middle].abs()) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= (256 * rank_id) - ranks_[rank_id].abs();

    uint64_t word_id = rank_id * 4;
    if (i < 128UL - ranks_[rank_id].rels[1]) {
      if (i >= 64UL - ranks_[rank_id].rels[0]) {
        word_id += 1;
        i -= 64UL - ranks_[rank_id].rels[0];
      }
    } else if (i < 192UL - ranks_[rank_id].rels[2]) {
      word_id += 2;
      i -= 128UL - ranks_[rank_id].rels[1];
    } else {
      word_id += 3;
      i -= 192UL - ranks_[rank_id].rels[2];
    }
    return (word_id * 64) + __builtin_ctzll(
      _pdep_u64(1UL << i, ~words_[word_id]));
  }

  uint64_t select1(uint64_t i) const {
    const uint64_t block_id = i / 256;
    uint64_t begin = select1s_[block_id];
    uint64_t end = select1s_[block_id + 1] + 1;
    if (begin + 10 >= end) {
      while (i >= ranks_[begin + 1].abs()) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < ranks_[middle].abs()) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= ranks_[rank_id].abs();

    uint64_t word_id = rank_id * 4;
    if (i < ranks_[rank_id].rels[1]) {
      if (i >= ranks_[rank_id].rels[0]) {
        word_id += 1;
        i -= ranks_[rank_id].rels[0];
      }
    } else if (i < ranks_[rank_id].rels[2]) {
      word_id += 2;
      i -= ranks_[rank_id].rels[1];
    } else {
      word_id += 3;
      i -= ranks_[rank_id].rels[2];
    }
    return (word_id * 64) + __builtin_ctzll(
      _pdep_u64(1UL << i, words_[word_id]));
  }

 private:
  vector<uint64_t> words_;
  uint64_t n_bits_;
  struct Rank {
    uint32_t abs_hi;
    uint8_t abs_lo;
    uint8_t rels[3];

    uint64_t abs() const {
      return ((uint64_t)abs_hi << 8) | abs_lo;
    }
    void set_abs(uint64_t abs) {
      abs_hi = (uint32_t)(abs >> 8);
      abs_lo = (uint8_t)abs;
    }
  };
  vector<Rank> ranks_;
  vector<uint64_t> select0s_;
  vector<uint64_t> select1s_;
  uint64_t n_zeros_;
  uint64_t n_ones_;

  void set_(uint64_t i, uint64_t bit) {
    if (bit) {
      words_[i / 64] |= (1UL << (i % 64));
    } else {
      words_[i / 64] &= ~(1UL << (i % 64));
    }
  }
};

struct Want {
  uint64_t id;
  uint64_t bit;
  uint64_t rank;
  uint64_t select;
};

constexpr uint64_t N_RANKS = 1UL << 28;

static void benchmark(const vector<uint64_t> &words,
  const vector<Want> &wants,
  const vector<uint64_t> &rank_queries,
  const vector<uint64_t> &select0_queries,
  const vector<uint64_t> &select1_queries) {
  uint64_t n_bits = words.size() * 64;

  printf("%-6luM", n_bits >> 20);

  BitVector bit_vector;
  high_resolution_clock::time_point begin = high_resolution_clock::now();
  for (uint64_t i = 0; i < n_bits; ++i) {
    bit_vector.add((words[i / 64] >> (i % 64)) & 1);
  }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  printf("| %7.3f", elapsed / n_bits);

  begin = high_resolution_clock::now();
  bit_vector.build();
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  printf("| %7.3f", elapsed / n_bits);

  for (uint64_t i = 0; i < wants.size(); ++i) {
    // cout << "i = " << i << ", id = " << wants[i].id << ", rank = " << wants[i].rank << endl;
    // if (bit_vector.rank(wants[i].id) != wants[i].rank) {
    //   cout << "rank = " << bit_vector.rank(wants[i].id) << ", want = " << wants[i].rank << endl;
    // }
    assert(bit_vector.rank(wants[i].id) == wants[i].rank);
    if (wants[i].bit) {
      assert(bit_vector.select1(wants[i].rank) == wants[i].id);
    } else {
      assert(bit_vector.select0(wants[i].id - wants[i].rank) == wants[i].id);
    }
  }

  begin = high_resolution_clock::now();
  uint64_t total = 0;
  for (uint64_t i = 0; i < rank_queries.size(); ++i) {
    total += bit_vector.rank(rank_queries[i]);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  assert(total != 0);
  printf("| %7.3f", elapsed / rank_queries.size());

  begin = high_resolution_clock::now();
  total = 0;
  for (uint64_t i = 0; i < select0_queries.size(); ++i) {
    total += bit_vector.select0(select0_queries[i]);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  assert(total != 0);
  printf("| %7.3f", elapsed / select0_queries.size());
  // printf("| %7.3f", 0.0);

  begin = high_resolution_clock::now();
  total = 0;
  for (uint64_t i = 0; i < select1_queries.size(); ++i) {
    total += bit_vector.select1(select1_queries[i]);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  assert(total != 0);
  printf("| %7.3f", elapsed / select1_queries.size());
  // printf("| %7.3f", 0.0);

  printf("\n");
}

int main() {
  random_device rnd_seed_gen;
  mt19937_64 rnd_engine(rnd_seed_gen());

  printf("n_bits |     add|   build|    rank| select0| select1\n");
  printf("------:|-------:|-------:|-------:|-------:|-------:\n");

  for (uint64_t n_bits = 1UL << 20; n_bits <= (1UL << 34); n_bits <<= 1) {
    vector<uint64_t> words(n_bits / 64);
    for (uint64_t i = 0; i < words.size(); ++i) {
      words[i] = rnd_engine();
    }

    vector<Want> wants(words.size());
    uint64_t n_ones = 0;
    for (uint64_t i = 0; i < words.size(); ++i) {
      wants[i].id = (i * 64) + (rnd_engine() % 64);
      wants[i].bit = (words[i] >> (wants[i].id % 64)) & 1;
      wants[i].rank = n_ones + __builtin_popcountll(
        words[i] & ((1UL << (wants[i].id % 64)) - 1));
      n_ones += __builtin_popcountll(words[i]);
    }
    uint64_t n_zeros = n_bits - n_ones;

    vector<uint64_t> rank_queries(1UL << 26);
    for (uint64_t i = 0; i < rank_queries.size(); ++i) {
      rank_queries[i] = rnd_engine() % n_bits;
    }
    vector<uint64_t> select0_queries(1UL << 26);
    for (uint64_t i = 0; i < select0_queries.size(); ++i) {
      select0_queries[i] = rnd_engine() % n_zeros;
    }
    vector<uint64_t> select1_queries(1UL << 26);
    for (uint64_t i = 0; i < select1_queries.size(); ++i) {
      select1_queries[i] = rnd_engine() % n_ones;
    }

    // cout << "#bits = " << n_bits << ", #ones = " << n_ones << endl;
    benchmark(words, wants, rank_queries, select0_queries, select1_queries);
  }
}
