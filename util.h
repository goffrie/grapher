#ifndef _UTIL_H_
#define _UTIL_H_

inline bool isIntegral(Number a) {
    const Number b = floor(a + 0.5);
    const Number epsilon = 1e-4;
    return (a < b + epsilon) && (b < a + epsilon);
}

inline double rnd_d(double num) {
    return (num >= 0.0) ? floor(num + 0.5) : ceil(num - 0.5);
}

inline qreal roundOneDigit(qreal n) {
    if (n < 0) return -roundOneDigit(-n);
    if (n == 0) return 0;
    const qreal l = std::log10(n);
    const qreal f = std::pow(10, l - std::floor(l));
    return rnd_d(f) * std::pow(10, std::floor(l));
}

/*QString textRoundOneDigit(qreal n) {
    if (n < 0) return QString("-") + textRoundOneDigit(-n);
    if (n == 0) return QString("0");
    const qreal l = std::log10(n);
    const qreal fl = std::floor(l);
    long il = (long)fl;
    int digit = qRound(std::pow(10, l - fl));
    if (digit > 9) {
        digit /= 10;
        ++il;
    }
    const QString f = QString::number(digit);
    if (il >= 0) {
        return f + QString(il, '0');
    } else {
        return QString("0.") + QString(-il-1, '0') + f;
    }
}*/

#endif
