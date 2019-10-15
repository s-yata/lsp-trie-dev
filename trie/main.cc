#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

struct Level {
  vector<bool> louds;
  vector<bool> outs;
  vector<uint8_t> labels;

  Level() : louds(), outs(), labels() {}
};

struct Trie {
  vector<Level> levels;
  string last_key;
  uint64_t n_keys;
  uint64_t n_nodes;
  uint64_t n_pat_nodes;

  Trie() : levels(2), last_key(), n_keys(0), n_nodes(1), n_pat_nodes(1) {
    levels[0].louds.push_back(true);
    levels[0].louds.push_back(false);
    levels[1].louds.push_back(false);
    levels[0].outs.push_back(false);
    levels[0].labels.push_back(' ');
  }

  void add(const std::string &key) {
    assert(n_keys == 0 || key > last_key);
    if (key.empty()) {
      levels[0].outs.back() = true;
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
        assert(!levels[i + 1].louds.back());
        if (i != 0 && !levels[i].outs.back() &&
          levels[i + 1].louds.size() >= 2 &&
          levels[i + 1].louds[levels[i + 1].louds.size() - 2] &&
          (levels[i + 1].louds.size() == 2 ||
          !levels[i + 1].louds[levels[i + 1].louds.size() - 3])) {
          // New branch, not root and not out.
          ++n_pat_nodes;
        }
        levels[i + 1].louds.back() = true;
        levels[i + 1].louds.push_back(false);
        levels[i + 1].outs.push_back(false);
        levels[i + 1].labels.push_back(key[i]);
        ++n_nodes;
        break;
      }
    }
    for (++i; i < key.length(); ++i) {
      levels[i + 1].louds.push_back(true);
      levels[i + 1].louds.push_back(false);
      levels[i + 1].outs.push_back(false);
      levels[i + 1].labels.push_back(key[i]);
      ++n_nodes;
    }
    levels[key.length() + 1].louds.push_back(false);
    levels[key.length()].outs.back() = true;
    last_key = key;
    ++n_keys;
    ++n_pat_nodes;
  }
};

int main() {
  ios_base::sync_with_stdio(false);
  Trie trie;
  string line;
  while (getline(cin, line)) {
    trie.add(line);
  }
  cout << "#keys = " << trie.n_keys << endl;
  cout << "#nodes = " << trie.n_nodes << endl;
  cout << "#pat_nodes = " << trie.n_pat_nodes << endl;
  cout << "#levels = " << trie.levels.size() << endl;
  uint64_t louds_size = 0;
  uint64_t outs_size = 0;
  uint64_t labels_size = 0;
  for (uint64_t i = 0; i < trie.levels.size(); ++i) {
    const Level &level = trie.levels[i];
    louds_size += level.louds.size();
    outs_size += level.outs.size();
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
  return 0;
}
