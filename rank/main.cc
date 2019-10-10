#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using namespace std;
using namespace std::chrono;

// class BitVector {
//  public:
//   virtual ~BitVector() {}
//
//   virtual uint64_t rank(uint64_t i) const = 0;
//
//   virtual string name() const = 0;
//   virtual uint64_t size() const = 0;
// };

class BitVectorBase {
 public:
  BitVectorBase(const vector<uint64_t> &values)
    : values_(values), ranks_abs_(), ranks_rel_() {
    uint64_t n_blks = values_.size() / 8;
    ranks_abs_.resize(n_blks);
    ranks_rel_.resize(n_blks);
    uint64_t count = 0;
    for (uint64_t blk_id = 0; blk_id < n_blks; ++blk_id) {
      ranks_abs_[blk_id] = count;
      for (uint64_t j = 0; j < 8; ++j) {
        ranks_rel_[blk_id] |= (count - ranks_abs_[blk_id]) << (63 - (9 * j));
        count += __builtin_popcountll(values_[(blk_id * 8) + j]);
      }
    }
  }
  ~BitVectorBase() {}

  uint64_t rank(uint64_t i) const {
    uint64_t div = i / 64;
    uint64_t mod = i % 64;
    uint64_t abs = div / 8;
    uint64_t rel = div % 8;
    uint64_t n = ranks_abs_[abs];
    n += (ranks_rel_[abs] >> (63 - (9 * rel))) & 0x1FFUL;
    n += __builtin_popcountll(values_[div] & ((1UL << mod) - 1));
    return n;
  }

  string name() const {
    return "Base";
  }
  uint64_t size() const {
    return sizeof(uint64_t) * values_.size()
      + sizeof(uint64_t) * ranks_abs_.size()
      + sizeof(uint64_t) * ranks_rel_.size();
  }

 private:
  const vector<uint64_t> values_;
  vector<uint64_t> ranks_abs_;
  vector<uint64_t> ranks_rel_;
};

class BitVectorPacked64 {
 public:
  BitVectorPacked64(const vector<uint64_t> &values)
    : values_(values), ranks_() {
    uint64_t n_blks = values_.size() / 8;
    ranks_.resize(n_blks);
    uint64_t count = 0;
    for (uint64_t blk_id = 0; blk_id < n_blks; ++blk_id) {
      ranks_[blk_id].abs = count;
      for (uint64_t j = 0; j < 8; ++j) {
        ranks_[blk_id].rel |= (count - ranks_[blk_id].abs) << (63 - (9 * j));
        count += __builtin_popcountll(values_[(blk_id * 8) + j]);
      }
    }
  }
  ~BitVectorPacked64() {}

  uint64_t rank(uint64_t i) const {
    uint64_t div = i / 64;
    uint64_t mod = i % 64;
    uint64_t abs = div / 8;
    uint64_t rel = div % 8;
    uint64_t n = ranks_[abs].abs;
    n += (ranks_[abs].rel >> (63 - (9 * rel))) & 0x1FFUL;
    n += __builtin_popcountll(values_[div] & ((1UL << mod) - 1));
    return n;
  }

  string name() const {
    return "Packed64";
  }
  uint64_t size() const {
    return sizeof(uint64_t) * values_.size()
      + sizeof(Rank) * ranks_.size();
  }

 private:
  struct Rank {
    uint64_t abs;
    uint64_t rel;
  };
  const vector<uint64_t> values_;
  vector<Rank> ranks_;
};

class BitVectorPacked32 {
 public:
  BitVectorPacked32(const vector<uint64_t> &values)
    : values_(values), ranks_(), ranks_abs_() {
    uint64_t n_blks = values_.size() / 8;
    ranks_.resize(n_blks);
    uint64_t n_abs_blks = (values_.size() >> 25) + 1;
    ranks_abs_.resize(n_abs_blks);
    uint64_t count = 0;
    uint64_t abs = 0;
    for (uint64_t blk_id = 0; blk_id < n_blks; ++blk_id) {
      if ((blk_id % (1UL << 22)) == 0) {
        abs = count;
        ranks_abs_[blk_id >> 22] = abs;
      }
      ranks_[blk_id].abs = count - abs;
      uint64_t rel = 0;
      for (uint64_t j = 0; j < 8; ++j) {
        rel |= (count - abs - ranks_[blk_id].abs) << (63 - (9 * j));
        if ((blk_id * 8) + j < values_.size()) {
          count += __builtin_popcountll(values_[(blk_id * 8) + j]);
        }
      }
      ranks_[blk_id].rel_hi = (uint32_t)(rel >> 32);
      ranks_[blk_id].rel_lo = (uint32_t)rel;
    }
  }
  ~BitVectorPacked32() {}

  uint64_t rank(uint64_t i) const {
    uint64_t div = i / 64;
    uint64_t mod = i % 64;
    uint64_t hi = div / 8;
    uint64_t lo = div % 8;
    uint64_t rel = ((uint64_t)ranks_[hi].rel_hi << 32) | ranks_[hi].rel_lo;
    uint64_t n = ranks_abs_[i >> 31];
    n += ranks_[hi].abs;
    n += (rel >> (63 - (9 * lo))) & 0x1FFUL;
    n += __builtin_popcountll(values_[div] & ((1UL << mod) - 1));
    return n;
  }

  string name() const {
    return "Packed32";
  }
  uint64_t size() const {
    return sizeof(uint64_t) * values_.size()
      + sizeof(Rank) * ranks_.size()
      + sizeof(uint64_t) * ranks_abs_.size();
  }

 private:
  struct Rank {
    uint32_t abs;
    uint32_t rel_lo;
    uint32_t rel_hi;
  };
  const vector<uint64_t> values_;
  vector<Rank> ranks_;
  vector<uint64_t> ranks_abs_;
};

class BitVectorBlock {
 public:
  BitVectorBlock(const vector<uint64_t> &values) : blocks_() {
    uint64_t n_blks = values.size() / 8;
    blocks_.resize(n_blks);
    uint64_t count = 0;
    for (uint64_t i = 0; i < n_blks; ++i) {
      blocks_[i].rank_abs = count;
      for (uint64_t j = 0; j < 8; ++j) {
        blocks_[i].values[j] = values[(i * 8) + j];
        blocks_[i].rank_rel |= (count - blocks_[i].rank_abs) << (63 - (9 * j));
        count += __builtin_popcountll(values[(i * 8) + j]);
      }
    }
  }
  ~BitVectorBlock() {}

  uint64_t rank(uint64_t i) const {
    uint64_t div = i / 64;
    uint64_t mod = i % 64;
    uint64_t lo = div % 8;
    uint64_t n = blocks_[i / 512].rank_abs;
    n += (blocks_[i / 512].rank_rel >> (63 - (9 * lo))) & 0x1FFUL;
    n += __builtin_popcountll(blocks_[i / 512].values[lo] & ((1UL << mod) - 1));
    return n;
  }

  string name() const {
    return "Block";
  }
  uint64_t size() const {
    return sizeof(Block) * blocks_.size();
  }

 private:
  struct Block {
    uint64_t values[8];
    uint64_t rank_abs;
    uint64_t rank_rel;
  };
  vector<Block> blocks_;
};

class BitVectorMarisa {
 public:
  BitVectorMarisa(const vector<uint64_t> &values)
    : values_(values), ranks_() {
    uint64_t n_blks = values_.size() / 8;
    ranks_.resize(n_blks);
    uint64_t count = 0;
    for (uint64_t blk_id = 0; blk_id < n_blks; ++blk_id) {
      ranks_[blk_id].abs = count;
      count += __builtin_popcountll(values_[(blk_id * 8)]);
      ranks_[blk_id].set_rel1(count - ranks_[blk_id].abs);
      count += __builtin_popcountll(values_[(blk_id * 8) + 1]);
      ranks_[blk_id].set_rel2(count - ranks_[blk_id].abs);
      count += __builtin_popcountll(values_[(blk_id * 8) + 2]);
      ranks_[blk_id].set_rel3(count - ranks_[blk_id].abs);
      count += __builtin_popcountll(values_[(blk_id * 8) + 3]);
      ranks_[blk_id].set_rel4(count - ranks_[blk_id].abs);
      count += __builtin_popcountll(values_[(blk_id * 8) + 4]);
      ranks_[blk_id].set_rel5(count - ranks_[blk_id].abs);
      count += __builtin_popcountll(values_[(blk_id * 8) + 5]);
      ranks_[blk_id].set_rel6(count - ranks_[blk_id].abs);
      count += __builtin_popcountll(values_[(blk_id * 8) + 6]);
      ranks_[blk_id].set_rel7(count - ranks_[blk_id].abs);
      count += __builtin_popcountll(values_[(blk_id * 8) + 7]);
    }
  }
  ~BitVectorMarisa() {}

  uint64_t rank(uint64_t i) const{
    const Rank &rank = ranks_[i / 512];
    uint64_t n = rank.abs;
    switch ((i / 64) % 8) {
      case 1: { n += rank.rel1(); break; }
      case 2: { n += rank.rel2(); break; }
      case 3: { n += rank.rel3(); break; }
      case 4: { n += rank.rel4(); break; }
      case 5: { n += rank.rel5(); break; }
      case 6: { n += rank.rel6(); break; }
      case 7: { n += rank.rel7(); break; }
    }
    n += __builtin_popcountll(values_[i / 64] & ((1UL << (i % 64)) - 1));
    return n;
  }

  string name() const {
    return "Marisa";
  }
  uint64_t size() const {
    return sizeof(uint64_t) * values_.size()
      + sizeof(Rank) * ranks_.size();
  }

 private:
  struct Rank {
    uint32_t abs;
    uint32_t rel_lo;
    uint32_t rel_hi;

    void set_rel1(uint64_t value) {
      rel_lo = (uint32_t)((rel_lo & ~0x7FU) | (value & 0x7FU));
    }
    void set_rel2(uint64_t value) {
      rel_lo = (uint32_t)((rel_lo & ~(0xFFU << 7)) | ((value & 0xFFU) << 7));
    }
    void set_rel3(uint64_t value) {
      rel_lo = (uint32_t)((rel_lo & ~(0xFFU << 15)) | ((value & 0xFFU) << 15));
    }
    void set_rel4(uint64_t value) {
      rel_lo = (uint32_t)((rel_lo & ~(0x1FFU << 23)) | ((value & 0x1FFU) << 23));
    }
    void set_rel5(uint64_t value) {
      rel_hi = (uint32_t)((rel_hi & ~0x1FFU) | (value & 0x1FFU));
    }
    void set_rel6(uint64_t value) {
      rel_hi = (uint32_t)((rel_hi & ~(0x1FFU << 9)) | ((value & 0x1FFU) << 9));
    }
    void set_rel7(uint64_t value) {
      rel_hi = (uint32_t)((rel_hi & ~(0x1FFU << 18)) | ((value & 0x1FFU) << 18));
    }

    uint64_t rel1() const {
      return rel_lo & 0x7FU;
    }
    uint64_t rel2() const {
      return (rel_lo >> 7) & 0xFFU;
    }
    uint64_t rel3() const {
      return (rel_lo >> 15) & 0xFFU;
    }
    uint64_t rel4() const {
      return (rel_lo >> 23) & 0x1FFU;
    }
    uint64_t rel5() const {
      return rel_hi & 0x1FFU;
    }
    uint64_t rel6() const {
      return (rel_hi >> 9) & 0x1FFU;
    }
    uint64_t rel7() const {
      return (rel_hi >> 18) & 0x1FFU;
    }
  };
  const vector<uint64_t> values_;
  vector<Rank> ranks_;
};

struct Pair {
  uint64_t id;
  uint64_t rank;
};

constexpr uint64_t N_RANKS = 1UL << 28;

template <typename T>
static void benchmark(const vector<uint64_t> &values,
  const vector<Pair> &pairs, const vector<uint64_t> &ids) {
  uint64_t n_bits = values.size() * 64;

  T bv(values);
  cout << " " << bv.name() << ":" << flush;
  for (uint64_t i = 0; i < pairs.size(); ++i) {
    assert(bv.rank(pairs[i].id) == pairs[i].rank);
  }
  cout << " size = " << bv.size() << flush;

  high_resolution_clock::time_point begin = high_resolution_clock::now();
  uint64_t total = 0, n = 0;
  while (n < N_RANKS) {
    uint64_t m = min(n_bits, N_RANKS - n);
    for (uint64_t i = 0; i < m; ++i) {
      total += bv.rank(i);
    }
    n += m;
  }
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  assert(total != 0);
  cout << ", seq. = " << (elapsed / n) << " ns" << flush;

  begin = high_resolution_clock::now();
  total = 0;
  for (uint64_t i = 0; i < ids.size(); ++i) {
    total += bv.rank(ids[i]);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  assert(total != 0);
  cout << ", rnd. = " << (elapsed / ids.size()) << " ns" << endl;
}

int main() {
  random_device rnd_seed_gen;
  mt19937_64 rnd_engine(rnd_seed_gen());

  for (uint64_t n_bits = 1UL << 18; n_bits < (1UL << 38); n_bits <<= 1) {
    vector<uint64_t> values(n_bits / 64);
    for (uint64_t i = 0; i < values.size(); ++i) {
      values[i] = rnd_engine();
    }
    vector<Pair> pairs(values.size());
    uint64_t n_ones = 0;
    for (uint64_t i = 0; i < values.size(); ++i) {
      pairs[i].id = (i * 64) + (rnd_engine() % 64);
      pairs[i].rank = n_ones + __builtin_popcountll(
        values[pairs[i].id / 64] & ((1UL << (pairs[i].id % 64)) - 1));
      n_ones += __builtin_popcountll(values[i]);
    }
    uint64_t n_ids = min(n_bits, 1UL << 26);
    vector<uint64_t> ids(n_ids);
    for (uint64_t i = 0; i < ids.size(); ++i) {
      ids[i] = rnd_engine() & (n_bits - 1);
    }

    cout << "#bits = " << n_bits << ", #ones = " << n_ones << endl;
    benchmark<BitVectorBase>(values, pairs, ids);
    benchmark<BitVectorPacked64>(values, pairs, ids);
    benchmark<BitVectorPacked32>(values, pairs, ids);
    benchmark<BitVectorBlock>(values, pairs, ids);
    if (n_bits < (1UL << 32)) {
      benchmark<BitVectorMarisa>(values, pairs, ids);
    }
  }
}
