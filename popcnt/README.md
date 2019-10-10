# popcnt

- CPU: Intel(R) Core(TM) i7-6800K CPU @ 3.40GHz

Function                                 |Speed
-----------------------------------------|-----------------
`__builtin_popcountll`                   |0.36967 values/ns
`__builtin_popcountll` + `-march=native` |2.16893 values/ns
`my_popcnt`                              |0.64247 values/ns

```cpp
inline uint64_t my_popcnt(uint64_t x) {
  x = (x & 0x5555555555555555ULL) + ((x & 0xAAAAAAAAAAAAAAAAULL) >> 1);
  x = (x & 0x3333333333333333ULL) + ((x & 0xCCCCCCCCCCCCCCCCULL) >> 2);
  x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x & 0xF0F0F0F0F0F0F0F0ULL) >> 4);
  x *= 0x0101010101010101ULL;
  return (x >> 56) & 0xFFULL;
}
```
