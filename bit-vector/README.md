# BitVector

- CPU: Intel(R) Core(TM) i7-6800K CPU @ 3.40GHz

## bit-vector-64 (nanoseconds)

n_bits |     add|   build|    rank| select0| select1
------:|-------:|-------:|-------:|-------:|-------:
1     M|   2.183|   0.097|   2.354|  20.603|  18.551
2     M|   2.005|   0.096|   4.109|  22.423|  20.342
4     M|   1.991|   0.094|   4.878|  24.952|  22.123
8     M|   1.985|   0.094|   5.206|  28.314|  25.555
16    M|   1.987|   0.093|   5.482|  30.824|  27.132
32    M|   1.985|   0.095|   5.676|  32.925|  28.711
64    M|   1.984|   0.096|   6.817|  34.629|  30.447
128   M|   1.988|   0.096|  14.415|  42.317|  36.010
256   M|   2.150|   0.112|  17.811|  54.913|  46.856
512   M|   2.249|   0.113|  21.259|  74.146|  62.587
1024  M|   2.289|   0.113|  23.513|  91.508|  78.424
2048  M|   2.306|   0.116|  25.061| 103.165|  92.610
4096  M|   2.326|   0.117|  25.953| 110.658|  99.470
8192  M|   2.340|   0.116|  26.604| 116.029| 104.943
16384 M|   2.333|   0.116|  29.249| 119.964| 108.592

## bit-vector-40 (nanoseconds)

n_bits |     add|   build|    rank| select0| select1
------:|-------:|-------:|-------:|-------:|-------:
1     M|   2.191|   0.151|   3.302|  16.720|  16.113
2     M|   2.022|   0.148|   4.426|  18.975|  17.366
4     M|   2.004|   0.149|   5.215|  21.251|  19.808
8     M|   1.998|   0.147|   5.890|  24.743|  22.642
16    M|   2.001|   0.148|   6.429|  26.549|  24.099
32    M|   2.001|   0.148|   6.790|  28.131|  24.934
64    M|   2.001|   0.149|   7.851|  29.470|  26.317
128   M|   2.003|   0.149|  14.712|  37.191|  33.611
256   M|   2.163|   0.176|  18.001|  52.788|  48.149
512   M|   2.257|   0.178|  22.031|  73.486|  64.588
1024  M|   2.295|   0.175|  25.011|  87.206|  75.867
2048  M|   2.366|   0.176|  26.878|  96.919|  82.486
4096  M|   2.375|   0.182|  28.163| 102.889|  86.713
8192  M|   2.393|   0.181|  28.811| 106.907|  88.961
16384 M|   2.349|   0.182|  32.828| 110.326|  92.337
