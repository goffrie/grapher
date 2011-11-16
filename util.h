#ifndef _UTIL_H_
#define _UTIL_H_

inline bool isIntegral(Number a) {
    const Number b = floor(a + 0.5);
    const Number epsilon = 1e-4;
    return (a < b + epsilon) && (b < a + epsilon);
}

inline int rnd(Number a) {
    return static_cast<int>(a + 0.5);
}

#endif