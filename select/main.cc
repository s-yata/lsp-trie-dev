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

// class BitVector {
//  public:
//   virtual ~BitVector() {}
//
//   virtual uint64_t select(uint64_t i) const = 0;
//
//   virtual string name() const = 0;
//   virtual uint64_t size() const = 0;
// };

const uint8_t SELECT_TABLE[8][256] = {
  {
    7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
  },
  {
    7, 7, 7, 1, 7, 2, 2, 1, 7, 3, 3, 1, 3, 2, 2, 1,
    7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1,
    5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1,
    6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1,
    5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    7, 7, 7, 1, 7, 2, 2, 1, 7, 3, 3, 1, 3, 2, 2, 1,
    7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    7, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1,
    5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    7, 6, 6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1,
    6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1,
    6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1,
    5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1
  },
  {
    7, 7, 7, 7, 7, 7, 7, 2, 7, 7, 7, 3, 7, 3, 3, 2,
    7, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2,
    7, 7, 7, 5, 7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2,
    7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,
    7, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2,
    7, 6, 6, 4, 6, 4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2,
    7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2,
    6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,
    7, 7, 7, 7, 7, 7, 7, 2, 7, 7, 7, 3, 7, 3, 3, 2,
    7, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2,
    7, 7, 7, 5, 7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2,
    7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2,
    7, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2,
    7, 6, 6, 4, 6, 4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2,
    7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2,
    6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2
  },
  {
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3,
    7, 7, 7, 7, 7, 7, 7, 4, 7, 7, 7, 4, 7, 4, 4, 3,
    7, 7, 7, 7, 7, 7, 7, 5, 7, 7, 7, 5, 7, 5, 5, 3,
    7, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3,
    7, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 6, 7, 6, 6, 3,
    7, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4, 6, 4, 4, 3,
    7, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3,
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 3,
    7, 7, 7, 7, 7, 7, 7, 4, 7, 7, 7, 4, 7, 4, 4, 3,
    7, 7, 7, 7, 7, 7, 7, 5, 7, 7, 7, 5, 7, 5, 5, 3,
    7, 7, 7, 5, 7, 5, 5, 4, 7, 5, 5, 4, 5, 4, 4, 3,
    7, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 6, 7, 6, 6, 3,
    7, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4, 6, 4, 4, 3,
    7, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 3,
    7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3
  },
  {
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 4,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5,
    7, 7, 7, 7, 7, 7, 7, 5, 7, 7, 7, 5, 7, 5, 5, 4,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6,
    7, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 6, 7, 6, 6, 4,
    7, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 6, 7, 6, 6, 5,
    7, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 4,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5,
    7, 7, 7, 7, 7, 7, 7, 5, 7, 7, 7, 5, 7, 5, 5, 4,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6,
    7, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 6, 7, 6, 6, 4,
    7, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 6, 7, 6, 6, 5,
    7, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6, 5, 6, 5, 5, 4
  },
  {
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6,
    7, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 6, 7, 6, 6, 5,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6,
    7, 7, 7, 7, 7, 7, 7, 6, 7, 7, 7, 6, 7, 6, 6, 5
  },
  {
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 6
  },
  {
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
  }
};

const uint64_t MASK_01 = 0x0101010101010101ULL;
const uint64_t MASK_0F = 0x0F0F0F0F0F0F0F0FULL;
const uint64_t MASK_33 = 0x3333333333333333ULL;
const uint64_t MASK_55 = 0x5555555555555555ULL;
const uint64_t MASK_80 = 0x8080808080808080ULL;

inline uint64_t select_bit(uint64_t i, uint64_t bit_id, uint64_t unit) {
  // return bit_id + __builtin_ctzll(_pdep_u64(1UL << i, unit));

  uint64_t counts;
  counts = unit - ((unit >> 1) & MASK_55);
  counts = (counts & MASK_33) + ((counts >> 2) & MASK_33);
  counts = (counts + (counts >> 4)) & MASK_0F;
  counts *= MASK_01;

  const uint64_t x = (counts | MASK_80) - ((i + 1) * MASK_01);
  const int skip = __builtin_ctzll((x & MASK_80) >> 7);

  bit_id += static_cast<uint64_t>(skip);
  unit >>= skip;
  i -= ((counts << 8) >> skip) & 0xFF;

  return bit_id + SELECT_TABLE[i][unit & 0xFF];
}

class BitVectorBase {
 public:
  BitVectorBase(const vector<uint64_t> &values)
    : values_(values), ranks_abs_(), ranks_rel_(), selects_() {
    uint64_t n_blks = values_.size() / 8;
    ranks_abs_.resize(n_blks + 1);
    ranks_rel_.resize(n_blks + 1);
    uint64_t n_ones = 0;
    for (uint64_t blk_id = 0; blk_id < n_blks; ++blk_id) {
      ranks_abs_[blk_id] = n_ones;
      for (uint64_t j = 0; j < 8; ++j) {
        ranks_rel_[blk_id] |= (n_ones - ranks_abs_[blk_id]) << (63 - (9 * j));
        uint64_t value_id = (blk_id * 8) + j;
        uint64_t value = values_[value_id];
        uint64_t n_pops = __builtin_popcountll(value);
        if (((n_ones - 1) / 512) != (((n_ones + n_pops) - 1) / 512)) {
          uint64_t count = n_ones;
          while (value != 0) {
            uint64_t pos = __builtin_ctzll(value);
            if (count % 512 == 0) {
              selects_.push_back((value_id * 64) + pos);
              break;
            }
            value ^= 1UL << pos;
            ++count;
          }
        }
        n_ones += n_pops;
      }
    }
    ranks_abs_.back() = n_ones;
    selects_.push_back(values_.size() * 64);
  }
  ~BitVectorBase() {}

  uint64_t select(uint64_t i) const {
    const uint64_t blk_id = i / 512;
    if ((i % 512) == 0) {
      return selects_[blk_id];
    }
    uint64_t begin = selects_[blk_id] / 512;
    uint64_t end = (selects_[blk_id + 1] + 511) / 512;
    if (begin + 10 >= end) {
      while (i >= ranks_abs_[begin + 1]) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < ranks_abs_[middle]) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= ranks_abs_[rank_id];

    uint64_t rel = ranks_rel_[rank_id];
    uint64_t value_id = rank_id * 8;
    if (i < ((rel >> (63 - (9 * 4))) & 0x1FFUL)) {
      if (i < ((rel >> (63 - (9 * 2))) & 0x1FFUL)) {
        if (i >= ((rel >> (63 - (9 * 1))) & 0x1FFUL)) {
          value_id += 1;
          i -= (rel >> (63 - (9 * 1))) & 0x1FFUL;
        }
      } else if (i < ((rel >> (63 - (9 * 3))) & 0x1FFUL)) {
        value_id += 2;
        i -= ((rel >> (63 - (9 * 2))) & 0x1FFUL);
      } else {
        value_id += 3;
        i -= ((rel >> (63 - (9 * 3))) & 0x1FFUL);
      }
    } else if (i < ((rel >> (63 - (9 * 6))) & 0x1FFUL)) {
      if (i < ((rel >> (63 - (9 * 5))) & 0x1FFUL)) {
        value_id += 4;
        i -= ((rel >> (63 - (9 * 4))) & 0x1FFUL);
      } else {
        value_id += 5;
        i -= ((rel >> (63 - (9 * 5))) & 0x1FFUL);
      }
    } else if (i < ((rel >> (63 - (9 * 7))) & 0x1FFUL)) {
      value_id += 6;
      i -= ((rel >> (63 - (9 * 6))) & 0x1FFUL);
    } else {
      value_id += 7;
      i -= ((rel >> (63 - (9 * 7))) & 0x1FFUL);
    }

    return select_bit(i, value_id * 64, values_[value_id]);
  }

  string name() const {
    return "Base";
  }
  uint64_t size() const {
    return sizeof(uint64_t) * values_.size()
      + sizeof(uint64_t) * ranks_abs_.size()
      + sizeof(uint64_t) * ranks_rel_.size()
      + sizeof(uint64_t) * selects_.size();
  }

 private:
  const vector<uint64_t> values_;
  vector<uint64_t> ranks_abs_;
  vector<uint64_t> ranks_rel_;
  vector<uint64_t> selects_;
};

class BitVectorPacked64 {
 public:
  BitVectorPacked64(const vector<uint64_t> &values)
    : values_(values), ranks_(), selects_() {
    uint64_t n_blks = values_.size() / 8;
    ranks_.resize(n_blks + 1);
    uint64_t n_ones = 0;
    for (uint64_t blk_id = 0; blk_id < n_blks; ++blk_id) {
      ranks_[blk_id].abs = n_ones;
      for (uint64_t j = 0; j < 8; ++j) {
        ranks_[blk_id].rel |= (n_ones - ranks_[blk_id].abs) << (63 - (9 * j));
        uint64_t value_id = (blk_id * 8) + j;
        uint64_t value = values_[value_id];
        uint64_t n_pops = __builtin_popcountll(value);
        if (((n_ones - 1) / 512) != (((n_ones + n_pops) - 1) / 512)) {
          uint64_t count = n_ones;
          while (value != 0) {
            uint64_t pos = __builtin_ctzll(value);
            if (count % 512 == 0) {
              selects_.push_back((value_id * 64) + pos);
              break;
            }
            value ^= 1UL << pos;
            ++count;
          }
        }
        n_ones += n_pops;
      }
    }
    ranks_.back().abs = n_ones;
    selects_.push_back(values_.size() * 64);
  }
  ~BitVectorPacked64() {}

  uint64_t select(uint64_t i) const {
    const uint64_t blk_id = i / 512;
    if ((i % 512) == 0) {
      return selects_[blk_id];
    }
    uint64_t begin = selects_[blk_id] / 512;
    uint64_t end = (selects_[blk_id + 1] + 511) / 512;
    if (begin + 10 >= end) {
      while (i >= ranks_[begin + 1].abs) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < ranks_[middle].abs) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= ranks_[rank_id].abs;

    uint64_t rel = ranks_[rank_id].rel;
    uint64_t value_id = rank_id * 8;
    if (i < ((rel >> (63 - (9 * 4))) & 0x1FFUL)) {
      if (i < ((rel >> (63 - (9 * 2))) & 0x1FFUL)) {
        if (i >= ((rel >> (63 - (9 * 1))) & 0x1FFUL)) {
          value_id += 1;
          i -= (rel >> (63 - (9 * 1))) & 0x1FFUL;
        }
      } else if (i < ((rel >> (63 - (9 * 3))) & 0x1FFUL)) {
        value_id += 2;
        i -= ((rel >> (63 - (9 * 2))) & 0x1FFUL);
      } else {
        value_id += 3;
        i -= ((rel >> (63 - (9 * 3))) & 0x1FFUL);
      }
    } else if (i < ((rel >> (63 - (9 * 6))) & 0x1FFUL)) {
      if (i < ((rel >> (63 - (9 * 5))) & 0x1FFUL)) {
        value_id += 4;
        i -= ((rel >> (63 - (9 * 4))) & 0x1FFUL);
      } else {
        value_id += 5;
        i -= ((rel >> (63 - (9 * 5))) & 0x1FFUL);
      }
    } else if (i < ((rel >> (63 - (9 * 7))) & 0x1FFUL)) {
      value_id += 6;
      i -= ((rel >> (63 - (9 * 6))) & 0x1FFUL);
    } else {
      value_id += 7;
      i -= ((rel >> (63 - (9 * 7))) & 0x1FFUL);
    }

    return select_bit(i, value_id * 64, values_[value_id]);
  }

  string name() const {
    return "Packed64";
  }
  uint64_t size() const {
    return sizeof(uint64_t) * values_.size()
      + sizeof(Rank) * ranks_.size()
      + sizeof(uint64_t) * selects_.size();
  }

 private:
  struct Rank {
    uint64_t abs;
    uint64_t rel;
  };
  const vector<uint64_t> values_;
  vector<Rank> ranks_;
  vector<uint64_t> selects_;
};

class BitVectorBlock {
 public:
  BitVectorBlock(const vector<uint64_t> &values) : blocks_(), selects_() {
    uint64_t n_blks = values.size() / 8;
    blocks_.resize(n_blks + 1);
    uint64_t n_ones = 0;
    for (uint64_t blk_id = 0; blk_id < n_blks; ++blk_id) {
      blocks_[blk_id].rank_abs = n_ones;
      for (uint64_t j = 0; j < 8; ++j) {
        uint64_t value_id = (blk_id * 8) + j;
        uint64_t value = values[value_id];
        blocks_[blk_id].values[j] = value;
        blocks_[blk_id].rank_rel |=
          (n_ones - blocks_[blk_id].rank_abs) << (63 - (9 * j));
        uint64_t n_pops = __builtin_popcountll(value);
        if (((n_ones - 1) / 512) != (((n_ones + n_pops) - 1) / 512)) {
          uint64_t count = n_ones;
          while (value != 0) {
            uint64_t pos = __builtin_ctzll(value);
            if (count % 512 == 0) {
              selects_.push_back((value_id * 64) + pos);
              break;
            }
            value ^= 1UL << pos;
            ++count;
          }
        }
        n_ones += n_pops;
      }
    }
    blocks_.back().rank_abs = n_ones;
    selects_.push_back(values.size() * 64);
  }
  ~BitVectorBlock() {}

  uint64_t select(uint64_t i) const {
    const uint64_t blk_id = i / 512;
    if ((i % 512) == 0) {
      return selects_[blk_id];
    }
    uint64_t begin = selects_[blk_id] / 512;
    uint64_t end = (selects_[blk_id + 1] + 511) / 512;
    if (begin + 10 >= end) {
      while (i >= blocks_[begin + 1].rank_abs) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < blocks_[middle].rank_abs) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= blocks_[rank_id].rank_abs;

    uint64_t rel = blocks_[rank_id].rank_rel;
    uint64_t value_id = rank_id * 8;
    if (i < ((rel >> (63 - (9 * 4))) & 0x1FFUL)) {
      if (i < ((rel >> (63 - (9 * 2))) & 0x1FFUL)) {
        if (i >= ((rel >> (63 - (9 * 1))) & 0x1FFUL)) {
          value_id += 1;
          i -= (rel >> (63 - (9 * 1))) & 0x1FFUL;
        }
      } else if (i < ((rel >> (63 - (9 * 3))) & 0x1FFUL)) {
        value_id += 2;
        i -= ((rel >> (63 - (9 * 2))) & 0x1FFUL);
      } else {
        value_id += 3;
        i -= ((rel >> (63 - (9 * 3))) & 0x1FFUL);
      }
    } else if (i < ((rel >> (63 - (9 * 6))) & 0x1FFUL)) {
      if (i < ((rel >> (63 - (9 * 5))) & 0x1FFUL)) {
        value_id += 4;
        i -= ((rel >> (63 - (9 * 4))) & 0x1FFUL);
      } else {
        value_id += 5;
        i -= ((rel >> (63 - (9 * 5))) & 0x1FFUL);
      }
    } else if (i < ((rel >> (63 - (9 * 7))) & 0x1FFUL)) {
      value_id += 6;
      i -= ((rel >> (63 - (9 * 6))) & 0x1FFUL);
    } else {
      value_id += 7;
      i -= ((rel >> (63 - (9 * 7))) & 0x1FFUL);
    }

    return select_bit(i, value_id * 64,
      blocks_[value_id / 8].values[value_id % 8]);
  }

  string name() const {
    return "Block";
  }
  uint64_t size() const {
    return sizeof(Block) * blocks_.size()
      + sizeof(uint64_t) * selects_.size();
  }

 private:
  struct Block {
    uint64_t values[8];
    uint64_t rank_abs;
    uint64_t rank_rel;
  };
  vector<Block> blocks_;
  vector<uint64_t> selects_;
};

struct Pair {
  uint64_t id;
  uint64_t select;
};

constexpr uint64_t N_SELECTS = 1UL << 24;

template <typename T>
static void benchmark(const vector<uint64_t> &values,
  const vector<Pair> &pairs, const vector<uint64_t> &ids, uint64_t n_ones) {
  uint64_t n_bits = values.size() * 64;

  high_resolution_clock::time_point begin = high_resolution_clock::now();
  T bv(values);
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  printf(" %s:", bv.name().c_str());
  // printf(" size = %lu,", bv.size());
  printf(" bld = %.3f ns,", elapsed / n_bits);

  for (uint64_t i = 0; i < pairs.size(); ++i) {
    assert(bv.select(pairs[i].id) == pairs[i].select);
  }

  begin = high_resolution_clock::now();
  uint64_t total = 0, n = 0;
  while (n < N_SELECTS) {
    uint64_t m = min(n_ones, N_SELECTS - n);
    for (uint64_t i = 0; i < m; ++i) {
      total += bv.select(i);
    }
    n += m;
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  assert(total != 0);
  printf(" seq = %.3f ns,", elapsed / n);

  begin = high_resolution_clock::now();
  total = 0;
  for (uint64_t i = 0; i < ids.size(); ++i) {
    total += bv.select(ids[i]);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  assert(total != 0);
  printf(" rnd = %.3f ns\n", elapsed / ids.size());
}

int main() {
  random_device rnd_seed_gen;
  mt19937_64 rnd_engine(rnd_seed_gen());

  for (uint64_t n_bits = 1UL << 18; n_bits < (1UL << 38); n_bits <<= 1) {
    vector<uint64_t> values(n_bits / 64);
    for (uint64_t i = 0; i < values.size(); ++i) {
      values[i] = rnd_engine();
    }

    vector<Pair> pairs;
    uint64_t n_ones = 0;
    for (uint64_t i = 0; i < values.size(); ++i) {
      uint64_t n = __builtin_popcountll(values[i]);
      if (n != 0) {
        uint64_t m = rnd_engine() % n;
        uint64_t value = values[i];
        for (uint64_t j = 0; j < m; ++j) {
          value ^= 1UL << __builtin_ctzll(value);
        }
        pairs.push_back(Pair{n_ones + m, (i * 64) + __builtin_ctzll(value)});
        n_ones += n;
      }
    }

    uint64_t n_ids = min(n_bits, 1UL << 24);
    vector<uint64_t> ids(n_ids);
    for (uint64_t i = 0; i < ids.size(); ++i) {
      ids[i] = rnd_engine() % n_ones;
    }

    cout << "#bits = " << n_bits << ", #ones = " << n_ones << endl;
    benchmark<BitVectorBase>(values, pairs, ids, n_ones);
    benchmark<BitVectorPacked64>(values, pairs, ids, n_ones);
    benchmark<BitVectorBlock>(values, pairs, ids, n_ones);
  }
}
