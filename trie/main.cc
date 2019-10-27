#include <x86intrin.h>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

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
    if (n_bits_ % 512 == 0) {
      words_.resize((n_bits_ + 512) / 64);
    }
    set_(n_bits_, bit);
    ++n_bits_;
  }

  uint64_t n_bits() const {
    return n_bits_;
  }
  uint64_t n_zeros() const {
    return n_zeros_;
  }
  uint64_t n_ones() const {
    return n_ones_;
  }

  void build() {
    uint64_t n_blocks = words_.size() / 8;
    ranks_.resize(n_blocks + 1);
    uint64_t n_zeros = 0;
    uint64_t n_ones = 0;
    for (uint64_t block_id = 0; block_id < n_blocks; ++block_id) {
      ranks_[block_id].abs = n_ones;
      for (uint64_t j = 0; j < 8; ++j) {
        uint64_t rel = n_ones - ranks_[block_id].abs;
        ranks_[block_id].rel |= rel << (63 - (9 * j));

        uint64_t word_id = (block_id * 8) + j;
        uint64_t n_pops = __builtin_popcountll(words_[word_id]);
        uint64_t new_n_zeros = n_zeros + 64 - n_pops;
        if (((n_zeros + 511) / 512) != ((new_n_zeros + 511) / 512)) {
          uint64_t count = n_zeros;
          uint64_t word = ~words_[word_id];
          while (word != 0) {
            uint64_t pos = __builtin_ctzll(word);
            if (count % 512 == 0) {
              select0s_.push_back((word_id * 64) + pos);
              break;
            }
            word ^= 1UL << pos;
            ++count;
          }
        }
        n_zeros = new_n_zeros;

        uint64_t new_n_ones = n_ones + n_pops;
        if (((n_ones + 511) / 512) != ((new_n_ones + 511) / 512)) {
          uint64_t count = n_ones;
          uint64_t word = words_[word_id];
          while (word != 0) {
            uint64_t pos = __builtin_ctzll(word);
            if (count % 512 == 0) {
              select1s_.push_back((word_id * 64) + pos);
              break;
            }
            word ^= 1UL << pos;
            ++count;
          }
        }
        n_ones = new_n_ones;
      }
    }
    ranks_.back().abs = n_ones;
    select1s_.push_back(words_.size() * 64);
    select0s_.push_back(words_.size() * 64);
  }

  uint64_t rank(uint64_t i) const {
    uint64_t word_id = i / 64;
    uint64_t bit_id = i % 64;
    uint64_t rank_id = word_id / 8;  // i / 512
    uint64_t rel_id = word_id % 8;
    uint64_t n = ranks_[rank_id].abs;
    n += (ranks_[rank_id].rel >> (63 - (9 * rel_id))) & 0x1FFUL;
    n += __builtin_popcountll(words_[word_id] & ((1UL << bit_id) - 1));
    return n;
  }

  uint64_t select0(uint64_t i) const {
    const uint64_t block_id = i / 512;
    if ((i % 512) == 0) {
      return select0s_[block_id];
    }
    uint64_t begin = select0s_[block_id] / 512;
    uint64_t end = (select0s_[block_id + 1] + 511) / 512;
    if (begin + 10 >= end) {
      while (i >= ((begin + 1) * 512) - ranks_[begin + 1].abs) {
        ++begin;
      }
    } else {
      while (begin + 1 < end) {
        const uint64_t middle = (begin + end) / 2;
        if (i < (middle * 512) - ranks_[middle].abs) {
          end = middle;
        } else {
          begin = middle;
        }
      }
    }
    const uint64_t rank_id = begin;
    i -= (rank_id * 512) - ranks_[rank_id].abs;

    uint64_t rel = ranks_[rank_id].rel;
    uint64_t word_id = rank_id * 8;
    if (i < (64 * 4) - ((rel >> (63 - (9 * 4))) & 0x1FFUL)) {
      if (i < (64 * 2) - ((rel >> (63 - (9 * 2))) & 0x1FFUL)) {
        if (i >= (64 * 1) - ((rel >> (63 - (9 * 1))) & 0x1FFUL)) {
          word_id += 1;
          i -= (64 * 1) - ((rel >> (63 - (9 * 1))) & 0x1FFUL);
        }
      } else if (i < (64 * 3) - ((rel >> (63 - (9 * 3))) & 0x1FFUL)) {
        word_id += 2;
        i -= (64 * 2) - ((rel >> (63 - (9 * 2))) & 0x1FFUL);
      } else {
        word_id += 3;
        i -= (64 * 3) - ((rel >> (63 - (9 * 3))) & 0x1FFUL);
      }
    } else if (i < (64 * 6) - ((rel >> (63 - (9 * 6))) & 0x1FFUL)) {
      if (i < (64 * 5) - ((rel >> (63 - (9 * 5))) & 0x1FFUL)) {
        word_id += 4;
        i -= (64 * 4) - ((rel >> (63 - (9 * 4))) & 0x1FFUL);
      } else {
        word_id += 5;
        i -= (64 * 5) - ((rel >> (63 - (9 * 5))) & 0x1FFUL);
      }
    } else if (i < (64 * 7) - ((rel >> (63 - (9 * 7))) & 0x1FFUL)) {
      word_id += 6;
      i -= (64 * 6) - ((rel >> (63 - (9 * 6))) & 0x1FFUL);
    } else {
      word_id += 7;
      i -= (64 * 7) - ((rel >> (63 - (9 * 7))) & 0x1FFUL);
    }
    return (word_id * 64) + __builtin_ctzll(
      _pdep_u64(1UL << i, ~words_[word_id]));
  }

  uint64_t select1(uint64_t i) const {
    const uint64_t block_id = i / 512;
    if ((i % 512) == 0) {
      return select1s_[block_id];
    }
    uint64_t begin = select1s_[block_id] / 512;
    uint64_t end = (select1s_[block_id + 1] + 511) / 512;
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
    uint64_t word_id = rank_id * 8;
    if (i < ((rel >> (63 - (9 * 4))) & 0x1FFUL)) {
      if (i < ((rel >> (63 - (9 * 2))) & 0x1FFUL)) {
        if (i >= ((rel >> (63 - (9 * 1))) & 0x1FFUL)) {
          word_id += 1;
          i -= (rel >> (63 - (9 * 1))) & 0x1FFUL;
        }
      } else if (i < ((rel >> (63 - (9 * 3))) & 0x1FFUL)) {
        word_id += 2;
        i -= ((rel >> (63 - (9 * 2))) & 0x1FFUL);
      } else {
        word_id += 3;
        i -= ((rel >> (63 - (9 * 3))) & 0x1FFUL);
      }
    } else if (i < ((rel >> (63 - (9 * 6))) & 0x1FFUL)) {
      if (i < ((rel >> (63 - (9 * 5))) & 0x1FFUL)) {
        word_id += 4;
        i -= ((rel >> (63 - (9 * 4))) & 0x1FFUL);
      } else {
        word_id += 5;
        i -= ((rel >> (63 - (9 * 5))) & 0x1FFUL);
      }
    } else if (i < ((rel >> (63 - (9 * 7))) & 0x1FFUL)) {
      word_id += 6;
      i -= ((rel >> (63 - (9 * 6))) & 0x1FFUL);
    } else {
      word_id += 7;
      i -= ((rel >> (63 - (9 * 7))) & 0x1FFUL);
    }
    return (word_id * 64) + __builtin_ctzll(
      _pdep_u64(1UL << i, words_[word_id]));
  }

 private:
  vector<uint64_t> words_;
  uint64_t n_bits_;
  struct Rank {
    uint64_t abs;
    uint64_t rel;
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

struct Level {
  BitVector louds;
  BitVector outs;
  vector<uint8_t> labels;
  uint64_t offset;

  Level() : louds(), outs(), labels(), offset(0) {}
};

struct Trie {
  vector<Level> levels;
  string last_key;
  uint64_t n_keys;
  uint64_t n_nodes;
  uint64_t n_pat_nodes;

  Trie() : levels(2), last_key(), n_keys(0), n_nodes(1), n_pat_nodes(1) {
    levels[0].louds.add(1);
    levels[0].louds.add(0);
    levels[1].louds.add(0);
    levels[0].outs.add(0);
    levels[0].labels.push_back(' ');
  }

  void add(const std::string &key) {
    assert(n_keys == 0 || key > last_key);
    if (key.empty()) {
      levels[0].outs.set(0, 1);
      ++n_keys;
      return;
    }
    if (key.length() + 1 >= levels.size()) {
      levels.resize(key.length() + 2);
    }
    uint64_t i = 0;
    for ( ; i < key.length(); ++i) {
      uint8_t byte = key[i];
      if ((i == last_key.length()) ||
          (byte != levels[i + 1].labels.back())) {
        assert(!levels[i + 1].louds.get(levels[i + 1].louds.n_bits() - 1));
        if (i != 0 && !levels[i].outs.get(levels[i].outs.n_bits() - 1) &&
          levels[i + 1].louds.n_bits() >= 2 &&
          levels[i + 1].louds.get(levels[i + 1].louds.n_bits() - 2) &&
          (levels[i + 1].louds.n_bits() == 2 ||
          !levels[i + 1].louds.get(levels[i + 1].louds.n_bits() - 3))) {
          // New branch, not root and not out.
          ++n_pat_nodes;
        }
        levels[i + 1].louds.set(levels[i + 1].louds.n_bits() - 1, 1);
        levels[i + 1].louds.add(0);
        levels[i + 1].outs.add(0);
        levels[i + 1].labels.push_back(key[i]);
        ++n_nodes;
        break;
      }
    }
    for (++i; i < key.length(); ++i) {
      levels[i + 1].louds.add(1);
      levels[i + 1].louds.add(0);
      levels[i + 1].outs.add(0);
      levels[i + 1].labels.push_back(key[i]);
      ++n_nodes;
    }
    levels[key.length() + 1].louds.add(0);
    levels[key.length()].outs.set(levels[key.length()].outs.n_bits() - 1, 1);
    last_key = key;
    ++n_keys;
    ++n_pat_nodes;
  }

  void build() {
    uint64_t offset = 0;
    for (uint64_t i = 0; i < levels.size(); ++i) {
      levels[i].louds.build();
      levels[i].offset = offset;
      offset += levels[i].outs.n_ones();
    }
  }

  int64_t lookup(const string &query) const {
    if (query.length() >= levels.size()) {
      return false;
    }
    uint64_t node_id = 0;
    uint64_t rank = 0;
    for (uint64_t i = 0; i < query.length(); ++i) {
      const Level &level = levels[i + 1];
      if (rank != 0) {
        node_id = level.louds.select0(rank - 1) + 1;
        rank = node_id - rank;
      } else {
        node_id = 0;
      }
      uint8_t byte = query[i];
      for ( ; ; ) {
        if (!level.louds.get(node_id) || level.labels[rank] > byte) {
          return -1;
        }
        if (level.labels[rank] == byte) {
          break;
        }
        ++node_id;
        ++rank;
      }
      // cout << i << ": node_id = " << node_id << ", rank = " << rank << endl;
    }
    const Level &level = levels[query.length()];
    if (!level.outs.get(rank)) {
      return false;
    }
    return level.offset + rank;
  }
};

int main() {
  ios_base::sync_with_stdio(false);
  vector<string> keys;
  string line;
  while (getline(cin, line)) {
    keys.push_back(line);
  }

  high_resolution_clock::time_point begin = high_resolution_clock::now();
  Trie trie;
  for (uint64_t i = 0; i < keys.size(); ++i) {
    trie.add(keys[i]);
  }
  trie.build();
  high_resolution_clock::time_point end = high_resolution_clock::now();
  double elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  printf("build = %7.3f ns\n", elapsed / keys.size());

  cout << "#keys = " << trie.n_keys << endl;
  cout << "#nodes = " << trie.n_nodes << endl;
  cout << "#pat_nodes = " << trie.n_pat_nodes << endl;
  cout << "#levels = " << trie.levels.size() << endl;
  uint64_t louds_size = 0;
  uint64_t outs_size = 0;
  uint64_t labels_size = 0;
  for (uint64_t i = 0; i < trie.levels.size(); ++i) {
    const Level &level = trie.levels[i];
    louds_size += level.louds.n_bits();
    outs_size += level.outs.n_bits();
    labels_size += level.labels.size();
    // cout << i << ": ";
    // cout << "louds = " << level.louds.size();
    // cout << ", outs = " << level.outs.size();
    // cout << ", labels = " << level.labels.size();
    // cout << endl;
  }
  cout << "louds = " << louds_size << " bits" << endl;
  cout << "outs = " << outs_size << " bits" << endl;
  cout << "labels = " << labels_size  << " bytes" << endl;
  // uint64_t n_bits = louds_size + outs_size + labels_size;
  // cout << "#bits = " << n_bits << ", #bytes = " << (n_bits / 8) << endl;

  begin = high_resolution_clock::now();
  for (uint64_t i = 0; i < keys.size(); ++i) {
    assert(trie.lookup(keys[i]) != -1);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  printf("seq. lookup = %7.3f ns\n", elapsed / keys.size());

  random_shuffle(keys.begin(), keys.end());

  begin = high_resolution_clock::now();
  for (uint64_t i = 0; i < keys.size(); ++i) {
    assert(trie.lookup(keys[i]) != -1);
  }
  end = high_resolution_clock::now();
  elapsed = (double)duration_cast<nanoseconds>(end - begin).count();
  printf("rnd. lookup = %7.3f ns\n", elapsed / keys.size());

  return 0;
}
