#include "Expression.h"

#include <bitset>

Vector Constant::evaluateVector(size_t size) const {
    VectorR r = new Number[size];
    Number _c = c;
    for (std::size_t i = 0; i < size; ++i) r[i] = _c;
    return r;
}

Vector External::evaluateVector(size_t size) const {
    VectorR r = new Number[size];
    for (std::size_t i = 0; i < size; ++i) r[i] = c[i];
    return r;
}

inline Number powi(Number a, int b) {
    Number ret = 1.;
    bool neg = b < 0;
    if (neg) b = -b;
    while (b > 0) {
        if (b & 1) ret *= a;
        b >>= 1;
        a *= a;
    }
    return neg ? 1./ret : ret;
}
Number PowInt::evaluate() const {
    return powi(a->evaluate(), b);
}
Number PowInt::evaluate(Number _a) const {
    return powi(_a, b);
}

Vector PowInt::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    int _b = b;
    VectorR ret = new Number[size];
    for (std::size_t i = 0; i < size; ++i) ret[i] = 1.;
    bool neg = _b < 0;
    if (neg) _b = -_b;
    while (_b > 0) {
        if (_b & 1) for (std::size_t i = 0; i < size; ++i) ret[i] *= _a[i];
        _b >>= 1;
        for (std::size_t i = 0; i < size; ++i) _a[i] *= _a[i];
    }
    if (neg) for (std::size_t i = 0; i < size; ++i) ret[i] = 1. / ret[i];
    delete[] _a;
    return ret;
}

Number PolyGamma::evaluate() const {
    return gsl_sf_psi_n(b, a->evaluate());
}

Number PolyGamma::evaluate(Number _a) const {
    return gsl_sf_psi_n(b, _a);
}

Vector PolyGamma::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    for (std::size_t i = 0; i < size; ++i) _a[i] = gsl_sf_psi_n(b, _a[i]);
    return _a;
}

Vector Pow::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    VectorR _b = b->evaluateVector(size);
    for (std::size_t i = 0; i < size; ++i) _a[i] = std::pow(_a[i], _b[i]);
    delete[] _b;
    return _a;
}

Vector Add::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    VectorR _b = b->evaluateVector(size);
    for (std::size_t i = 0; i < size; ++i) _a[i] += _b[i];
    delete[] _b;
    return _a;
}

Vector Sub::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    VectorR _b = b->evaluateVector(size);
    for (std::size_t i = 0; i < size; ++i) _a[i] -= _b[i];
    delete[] _b;
    return _a;
}

Vector Mul::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    VectorR _b = b->evaluateVector(size);
    for (std::size_t i = 0; i < size; ++i) _a[i] *= _b[i];
    delete[] _b;
    return _a;
}

Vector Div::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    VectorR _b = b->evaluateVector(size);
    for (std::size_t i = 0; i < size; ++i) _a[i] /= _b[i];
    delete[] _b;
    return _a;
}

Vector Neg::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    for (std::size_t i = 0; i < size; ++i) _a[i] = -_a[i];
    return _a;
}

Expression* UnaryOp::simplify() const {
    Expression* as = a->simplify();
    Constant* ac = dynamic_cast<Constant*>(as);
    if (ac != NULL) {
        Number acc = ac->c;
        delete ac;
        return new Constant(evaluate(acc));
    }
    return construct(as);
}

Expression* Neg::simplify() const {
    Expression* as = a->simplify();
    Constant* ac = dynamic_cast<Constant*>(as);
    if (ac != NULL) {
        Number acc = ac->c;
        delete ac;
        return new Constant(evaluate(acc));
    }
    Neg* an = dynamic_cast<Neg*>(ac);
    if (an != NULL) {
        Expression* r = an->a->copy();
        delete an;
        return r;
    }
    return new Neg(as);
}

Expression* Add::simplify() const {
    Expression* as = a->simplify();
    Expression* bs = b->simplify();
    Constant* ac = dynamic_cast<Constant*>(as);
    Constant* bc = dynamic_cast<Constant*>(bs);
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    if (ac && bc) {
        delete as;
        delete bs;
        return new Constant(acc + bcc);
    }
    // 0 + b = b
    if (ac && acc == 0) {
        delete as;
        return bs;
    }
    // a + 0 = a
    if (bc && bcc == 0) {
        delete bs;
        return as;
    }
    Neg* an = dynamic_cast<Neg*>(as);
    Neg* bn = dynamic_cast<Neg*>(bs);
    // (-a) + (-b) = -(a+b)
    if (an && bn) {
        Neg* r = new Neg(new Add(an->a->copy(), bn->a->copy()));
        delete as;
        delete bs;
        return r;
    }
    // (-a) + b = b-a
    if (an) {
        Sub* r = new Sub(bs, an->a->copy());
        delete as;
        return r;
    }
    // a + (-b) = a-b
    if (bn) {
        Sub* r = new Sub(as, bn->a->copy());
        delete bs;
        return r;
    }
    return new Add(as, bs);
}

Expression* Sub::simplify() const {
    Add* _r = new Add(a->copy(), new Neg(b->copy()));
    Expression* r = _r->simplify();
    delete _r;
    return r;
}

Expression* Mul::simplify() const {
    Expression* as = a->simplify();
    Expression* bs = b->simplify();
    Constant* ac = dynamic_cast<Constant*>(as);
    Constant* bc = dynamic_cast<Constant*>(bs);
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    // constant folding
    if (ac && bc) {
        delete as;
        delete bs;
        return new Constant(acc * bcc);
    }
    // 0 * b = 0
    // a * 0 = 0
    if ((ac && acc == 0) || (bc && bcc == 0)) {
        delete as;
        delete bs;
        return new Constant(0);
    }
    // 1 * b = b
    if (ac && acc == 1) {
        delete as;
        return bs;
    }
    // a * 1 = a
    if (bc && bcc == 1) {
        delete bs;
        return as;
    }
    Neg* an = dynamic_cast<Neg*>(as);
    Neg* bn = dynamic_cast<Neg*>(bs);
    // (-a) * (-b) = a * b
    if (an && bn) {
        Mul* _r = new Mul(an->a->copy(), bn->a->copy());
        delete as;
        delete bs;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    // (-a) * b = -(a * b)
    if (an) {
        Neg* _r = new Neg(new Mul(an->a->copy(), bs));
        delete as;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    // a * (-b) = -(a * b)
    if (an) {
        Neg* _r = new Neg(new Mul(as, bn->a->copy()));
        delete bs;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    Div* ad = dynamic_cast<Div*>(as);
    Div* bd = dynamic_cast<Div*>(bs);
    // (aa/ab) * (ba/bb) = (aa * ba) / (ab * bb)
    if (ad && bd) {
        Div* _r = new Div(new Mul(ad->a->copy(), bd->a->copy()), new Mul(ad->b->copy(), bd->b->copy()));
        delete as;
        delete bs;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    // (aa/ab) * b = (aa * b) / (ab)
    if (ad) {
        Div* _r = new Div(new Mul(ad->a->copy(), bs), ad->b->copy());
        delete as;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    // a * (ba/bb) = (a * ba) / (bb)
    if (bd) {
        Div* _r = new Div(new Mul(as, bd->a->copy()), bd->b->copy());
        delete bs;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }  
    return new Mul(as, bs);
}

Expression* Div::simplify() const {
    Expression* as = a->simplify();
    Expression* bs = b->simplify();
    Constant* ac = dynamic_cast<Constant*>(as);
    Constant* bc = dynamic_cast<Constant*>(bs);
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    // constant folding
    if (ac && bc) {
        delete as;
        delete bs;
        return new Constant(acc / bcc);
    }
    // 0 / b = 0
    if (ac && acc == 0) {
        delete as;
        delete bs;
        return new Constant(0);
    }
    // a / 0 = NaN
    if (bc && bcc == 0) {
        delete as;
        delete bs;
        return new Constant(std::numeric_limits<Number>::quiet_NaN());
    }
    // a / 1 = a
    if (bc && bcc == 1) {
        delete bs;
        return as;
    }
    Neg* an = dynamic_cast<Neg*>(as);
    Neg* bn = dynamic_cast<Neg*>(bs);
    // (-a) / (-b) = a / b
    if (an && bn) {
        Div* _r = new Div(an->a->copy(), bn->a->copy());
        delete as;
        delete bs;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    // (-a) / b = -(a * b)
    if (an) {
        Neg* _r = new Neg(new Div(an->a->copy(), bs));
        delete as;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    // a * (-b) = -(a * b)
    if (an) {
        Neg* _r = new Neg(new Div(as, bn->a->copy()));
        delete bs;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    Div* ad = dynamic_cast<Div*>(as);
    Div* bd = dynamic_cast<Div*>(bs);
    // (aa/ab) / (ba/bb) = (aa * bb) / (ab * ba)
    if (ad && bd) {
        Div* _r = new Div(new Mul(ad->a->copy(), bd->b->copy()), new Mul(ad->b->copy(), bd->a->copy()));
        delete as;
        delete bs;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    // (aa/ab) / b = aa / (ab * b)
    if (ad) {
        Div* _r = new Div(ad->a->copy(), new Mul(ad->b->copy(), bs));
        delete as;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }
    // a / (ba/bb) = (a * bb) / (ba)
    if (bd) {
        Div* _r = new Div(new Mul(as, bd->b->copy()), bd->a->copy());
        delete bs;
        Expression* r = _r->simplify();
        delete _r;
        return r;
    }  
    return new Div(as, bs);
}

Expression* Pow::simplify() const {
    Expression* as = a->simplify();
    Expression* bs = b->simplify();
    Constant* ac = dynamic_cast<Constant*>(as);
    Constant* bc = dynamic_cast<Constant*>(bs);
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    // constant folding
    if (ac && bc) {
        delete as;
        delete bs;
        return new Constant(std::pow(acc, bcc));
    }
    // a ^ 0 = 1; (0^0 = 1)
    // 1 ^ b = 1
    if ((bc && bcc == 0) || (ac && acc == 1)) {
        delete as;
        delete bs;
        return new Constant(1);
    }
    // 0 ^ b = 0
    if (ac && acc == 0) {
        delete as;
        delete bs;
        return new Constant(0);
    }
    // a ^ 1 = a
    if (bc && bcc == 1) {
        delete bs;
        return as;
    }
    return new Pow(as, bs);
}

Expression* PowInt::simplify() const {
    Expression* as = a->simplify();
    Constant* ac = dynamic_cast<Constant*>(as);
    Number acc = ac ? ac->c : 0;
    // constant folding
    if (ac) {
        delete as;
        return new Constant(powi(acc, b));
    }
    // a ^ 0 = 1; (0^0 = 1)
    // 1 ^ b = 1
    if ((b == 0) || (ac && acc == 1)) {
        delete as;
        return new Constant(1);
    }
    // 0 ^ b = 0
    if (ac && acc == 0) {
        delete as;
        return new Constant(0);
    }
    // a ^ 1 = a
    if (b == 1) {
        return as;
    }
    return new PowInt(as, b);
}

// de^u/dx = de^u/du * du/dx = e^u * du/dx
Expression* Exp::derivative(const Variable& var) const {
    return new Mul(copy(), a->derivative(var));
}
// dln(u)/dx = dln(u)/du * du/dx = 1/u * du/dx
Expression* Ln::derivative(const Variable& var) const {
    return new Mul(new Div(new Constant(1), a->copy()), a->derivative(var));
}
// dsqrt(u)/dx = dsqrt(u)/du * du/dx = (1/2*sqrt(u)) * du/dx = (du/dx) / (2 * sqrt(u))
Expression* Sqrt::derivative(const Variable& var) const {
    return new Div(a->derivative(var), new Mul(new Constant(2), new Sqrt(a->copy())));
}
// dsin(u)/dx = dsin(u)/du * du/dx = cos(u) * du/dx
Expression* Sin::derivative(const Variable& var) const {
    return new Mul(new Cos(a->copy()), a->derivative(var));
}
// dcos(u)/dx = dcos(u)/du * du/dx = -sin(u) * du/dx
Expression* Cos::derivative(const Variable& var) const {
    return new Neg(new Mul(new Sin(a->copy()), a->derivative(var)));
}
// dtan(u)/dx = dtan(u)/du * du/dx = (1 + tan^2(u)) * du/dx
Expression* Tan::derivative(const Variable& var) const {
    return new Mul(new Add(new Constant(1), new PowInt(new Tan(a->copy()), 2)), a->derivative(var));
}
// d asin(u)/dx = d asin(u) / du * du/dx = (du/dx)/sqrt(1 - u^2)
Expression* Asin::derivative(const Variable& var) const {
    return new Div(a->derivative(var), new Sqrt(new Sub(new Constant(1), new PowInt(a->copy(), 2))));
}
// d acos(u)/dx = d acos(u) / du * du/dx = -(du/dx)/sqrt(1 - u^2)
Expression* Acos::derivative(const Variable& var) const {
    return new Neg(new Div(a->derivative(var), new Sqrt(new Sub(new Constant(1), new PowInt(a->copy(), 2)))));
}
// d atan(u)/dx = d atan(u) / du * du/dx = (du/dx)/(1 + u^2)
Expression* Atan::derivative(const Variable& var) const {
    return new Div(a->derivative(var), new Add(new Constant(1), new PowInt(a->copy(), 2)));
}

// d(u+v)/dx = du/dx + dv/dx
Expression* Add::derivative(const Variable& var) const {
    return new Add(a->derivative(var), b->derivative(var));
}
// d(u-v)/dx = du/dx - dv/dx
Expression* Sub::derivative(const Variable& var) const {
    return new Sub(a->derivative(var), b->derivative(var));
}
// d(u*v)/dx = udv/dx + vdu/dx
Expression* Mul::derivative(const Variable& var) const {
    return new Add(new Mul(a->copy(), b->derivative(var)), new Mul(b->copy(), a->derivative(var)));
}
// d(u/v)/dx = (vdu/dx - udv/dx) / (u^2)
Expression* Div::derivative(const Variable& var) const {
    return new Div(
                new Sub(new Mul(b->copy(), a->derivative(var)), new Mul(a->copy(), b->derivative(var))),
                new PowInt(b->copy(), 2)
            );
}
// d(u^v)/dx = d/dx e^(vlogu) = d(e^(vlogu))/d(vlogu) d(vlogu)/dx
// = e^vlogu * (logu dv/dx + v dlogu/dx) = u^v * (logu dv/dx + v dlogu/du du/dx)
// = u^v * (logu dv/dx + v/u du/dx)
Expression* Pow::derivative(const Variable& var) const {
    return new Mul(
                copy(),
                new Add(
                    new Mul(new Ln(a->copy()), b->derivative(var)),
                    new Mul(new Div(b->copy(), a->copy()), a->derivative(var))
                )
            );
}
// d(u^a)/dx = d(u^b)/du * du/dx = bu^(b-1) * du/dx
Expression* PowInt::derivative(const Variable& var) const {
    return new Mul(
                new Mul(new Constant(b), new PowInt(a->copy(), b-1)),
                a->derivative(var)
            );
}

// dgamma(u)/dx = dgamma(u)/du * du/dx = gamma(u)*psi_0(u) * du/dx
Expression* Gamma::derivative(const Variable& var) const {
    return new Mul(new Mul(new Gamma(a->copy()), new PolyGamma(a->copy(), 0)), a->derivative(var));
}

// dpsi_b(u)/dx = dpsi_b(u)/du * du/dx = psi_(b+1)(u) * du/dx
Expression* PolyGamma::derivative(const Variable& var) const {
    return new Mul(new PolyGamma(a->copy(), b+1), a->derivative(var));
}
