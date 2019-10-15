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

  Trie() : levels(2), last_key(), n_keys(0) {
    levels[0].louds.push_back(true);
    levels[0].louds.push_back(false);
    levels[1].louds.push_back(false);
    levels[0].outs.push_back(false);
    levels[0].labels.push_back(' ');
  }

  void add(const std::string &key) {
    assert(n_keys == 0 || key > last_key);
    if (key.length() + 1 >= levels.size()) {
      levels.resize(key.length() + 2);
    }
    uint64_t i = 0;
    for ( ; i < key.length(); ++i) {
      uint8_t byte = key[i];
      if (levels[i + 1].labels.empty() ||
          (byte != levels[i + 1].labels.back())) {
        levels[i + 1].louds.back() = true;
        levels[i + 1].louds.push_back(false);
        levels[i + 1].outs.push_back(false);
        levels[i + 1].labels.push_back(key[i]);
        break;
      }
    }
    for (++i; i < key.length(); ++i) {
      levels[i + 1].louds.push_back(true);
      levels[i + 1].louds.push_back(false);
      levels[i + 1].outs.push_back(false);
      levels[i + 1].labels.push_back(key[i]);
    }
    levels[key.length() + 1].louds.push_back(false);
    levels[key.length()].outs.back() = true;
    last_key = key;
    ++n_keys;
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
  cout << "#levels = " << trie.levels.size() << endl;
  uint64_t n_bits = 0;
  for (uint64_t i = 0; i < 10 && i < trie.levels.size(); ++i) {
    cout << i << ": ";
    const Level &level = trie.levels[i];
    cout << "louds = " << level.louds.size();
    cout << ", outs = " << level.outs.size();
    cout << ", labels = " << level.labels.size();
    n_bits += level.louds.size();
    n_bits += level.outs.size();
    n_bits += level.labels.size() * 8;
    cout << endl;
  }
  cout << "#bits = " << n_bits << ", #bytes = " << (n_bits / 8) << endl;
  return 0;
}