#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct KeyInfo {
  uint64_t offset:40;
  uint64_t length:24;
};

class KeySet {
 public:
  KeySet()
    : max_n_keys_(1UL << 30), max_n_bytes_(1UL << 34), bytes_(), infos_() {}

  bool add(const string &key) {
    if (n_keys() >= max_n_keys_ || n_bytes() + key.size() > max_n_bytes_) {
      return false;
    }
    infos_.push_back(KeyInfo{ bytes_.size(), key.size() });
    for (uint64_t i = 0; i < key.size(); ++i) {
      bytes_.push_back(key[i]);
    }
    return true;
  }
  uint64_t n_keys() const {
    return infos_.size();
  }
  uint64_t n_bytes() const {
    return bytes_.size();
  }
  uint64_t size() const {
    return bytes_.size() + sizeof(uint64_t) * infos_.size();
  }
  void sort() {
    sort(0, n_keys(), 0);
  }
  void print() const {
    for (uint64_t i = 0; i < n_keys(); ++i) {
      cout.write((const char *)&bytes_[infos_[i].offset], infos_[i].length) << '\n';
    }
  }

 private:
  uint64_t max_n_keys_;
  uint64_t max_n_bytes_;
  vector<uint8_t> bytes_;
  vector<KeyInfo> infos_;

  enum {
    THRESHOLD = 10
  };

  int get_label(uint64_t i, uint64_t depth) const {
    return (depth < infos_[i].length) ? bytes_[infos_[i].offset + depth] : -1;
  }
  int median(uint64_t a, uint64_t b, uint64_t c, uint64_t depth) {
    const int x = get_label(a, depth);
    const int y = get_label(b, depth);
    const int z = get_label(c, depth);
    if (x < y) {
      if (y < z) {
        return y;
      } else if (x < z) {
        return z;
      }
      return x;
    } else if (x < z) {
      return x;
    } else if (y < z) {
      return z;
    }
    return y;
  }
  int compare(KeyInfo lhs, KeyInfo rhs, uint64_t depth) {
    for (uint64_t i = depth; i < lhs.length; ++i) {
      if (i == rhs.length) {
        return 1;
      }
      if (bytes_[lhs.offset + i] != bytes_[rhs.offset + i]) {
        return bytes_[lhs.offset + i] - bytes_[rhs.offset + i];
      }
    }
    return (int)(lhs.length - rhs.length);
  }
  uint64_t insertion_sort(uint64_t l, uint64_t r, uint64_t depth) {
    assert(l <= r);

    uint64_t count = 1;
    for (uint64_t i = l + 1; i < r; ++i) {
      int result = 0;
      for (uint64_t j = i; j > l; --j) {
        result = compare(infos_[j - 1], infos_[j], depth);
        if (result <= 0) {
          break;
        }
        swap(infos_[j - 1], infos_[j]);
      }
      if (result != 0) {
        ++count;
      }
    }
    return count;
  }
  uint64_t sort(uint64_t l, uint64_t r, uint64_t depth) {
    assert(l <= r);

    std::size_t count = 0;
    while ((r - l) > THRESHOLD) {
      uint64_t pl = l;
      uint64_t pr = r;
      uint64_t pivot_l = l;
      uint64_t pivot_r = r;

      const int pivot = median(l, (l + (r - l) / 2), (r - 1), depth);
      for ( ; ; ) {
        while (pl < pr) {
          const int label = get_label(pl, depth);
          if (label > pivot) {
            break;
          } else if (label == pivot) {
            swap(infos_[pl], infos_[pivot_l]);
            ++pivot_l;
          }
          ++pl;
        }
        while (pl < pr) {
          const int label = get_label(--pr, depth);
          if (label < pivot) {
            break;
          } else if (label == pivot) {
            swap(infos_[pr], infos_[--pivot_r]);
          }
        }
        if (pl >= pr) {
          break;
        }
        swap(infos_[pl], infos_[pr]);
        ++pl;
      }
      while (pivot_l > l) {
        swap(infos_[--pivot_l], infos_[--pl]);
      }
      while (pivot_r < r) {
        swap(infos_[pivot_r], infos_[pr]);
        ++pivot_r;
        ++pr;
      }

      if (((pl - l) > (pr - pl)) || ((r - pr) > (pr - pl))) {
        if ((pr - pl) == 1) {
          ++count;
        } else if ((pr - pl) > 1) {
          if (pivot == -1) {
            ++count;
          } else {
            count += sort(pl, pr, depth + 1);
          }
        }

        if ((pl - l) < (r - pr)) {
          if ((pl - l) == 1) {
            ++count;
          } else if ((pl - l) > 1) {
            count += sort(l, pl, depth);
          }
          l = pr;
        } else {
          if ((r - pr) == 1) {
            ++count;
          } else if ((r - pr) > 1) {
            count += sort(pr, r, depth);
          }
          r = pl;
        }
      } else {
        if ((pl - l) == 1) {
          ++count;
        } else if ((pl - l) > 1) {
          count += sort(l, pl, depth);
        }

        if ((r - pr) == 1) {
          ++count;
        } else if ((r - pr) > 1) {
          count += sort(pr, r, depth);
        }

        l = pl, r = pr;
        if ((pr - pl) == 1) {
          ++count;
        } else if ((pr - pl) > 1) {
          if (pivot == -1) {
            l = r;
            ++count;
          } else {
            ++depth;
          }
        }
      }
    }

    if ((r - l) > 1) {
      count += insertion_sort(l, r, depth);
    }
    return count;
  }
};

int main() {
  ios_base::sync_with_stdio(false);

  KeySet keyset;
  string line;
  while (getline(cin, line)) {
    if (!keyset.add(line)) {
      break;
    }
  }

  cout << "n_keys = " << keyset.n_keys() << endl;
  cout << "n_bytes = " << keyset.n_bytes() << endl;
  cout << "size = " << keyset.size() << endl;

  keyset.sort();
  // keyset.print();

  return 0;
}
