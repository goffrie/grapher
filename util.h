#ifndef _UTIL_H_
#define _UTIL_H_

inline constexpr bool isIntegral(Number a) {
    constexpr Number b = floor(a + 0.5);
    constexpr Number epsilon = 1e-4;
    return (a < b + epsilon) && (b < a + epsilon);
}

inline constexpr int rnd(Number a) {
    return static_cast<int>(a + 0.5);
}

#endif