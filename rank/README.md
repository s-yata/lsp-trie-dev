# BitVector.rank

- CPU: Intel(R) Core(TM) i7-6800K CPU @ 3.40GHz

## Random access (nanoseconds)

#bits          |Base   |Packed64 |Packed32 |Block  |Marisa
--------------:|------:|--------:|--------:|------:|------:
262,144        |  1.725|    1.972|    2.228|  2.076|  8.126
524,288        |  1.788|    2.007|    2.292|  2.208|  8.138
1,048,576      |  2.117|    2.336|    2.626|  2.547|  8.133
2,097,152      |  4.076|    4.147|    4.275|  5.664|  8.305
4,194,304      |  5.462|    4.792|    5.125|  5.078|  8.566
8,388,608      |  6.704|    5.161|    5.569|  5.412|  8.611
16,777,216     |  7.636|    5.465|    5.930|  5.570|  8.641
33,554,432     |  8.237|    5.657|    6.181|  5.674|  8.662
67,108,864     |  9.518|    7.140|    7.500|  7.030|  9.248
134,217,728    | 15.554|   14.533|   15.909| 15.129| 13.940
268,435,456    | 18.910|   17.917|   19.873| 18.752| 16.842
536,870,912    | 25.138|   21.002|   22.751| 20.713| 19.961
1,073,741,824  | 29.937|   23.556|   25.371| 21.903| 23.014
2,147,483,648  | 32.948|   25.102|   27.272| 22.630| 24.866
4,294,967,296  | 35.094|   26.003|   28.495| 22.903|    N/A
8,589,934,592  | 37.903|   26.638|   29.456| 23.566|    N/A
17,179,869,184 | 43.167|   20.606|   35.531| 27.135|    N/A
34,359,738,368 | 52.695|   39.174|   41.583| 33.095|    N/A
68,719,476,736 | 67.554|   51.647|   50.472| 37.780|    N/A
137,438,953,472| 84.020|   66.137|   62.419| 40.890|    N/A