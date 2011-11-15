#include "Expression.h"

#include <bitset>

#define UNARY_FUNCTION(Name, func, sfunc) \
Number Name::evaluate(Number _a) const { return sfunc(_a); } \
Number Name::evaluate() const { return sfunc(a->evaluate()); } \
Vector Name::evaluateVector(std::size_t size) const { \
    VectorR _a = a->evaluateVector(size); \
    for (std::size_t i = 0; i < size; ++i) _a[i] = sfunc(_a[i]); \
    return _a; \
}

UNARY_FUNCTION(Neg, -, -)
UNARY_FUNCTION(Exp, exp, std::exp)
UNARY_FUNCTION(Ln, log, std::log)
UNARY_FUNCTION(Sqrt, sqrt, std::sqrt)
UNARY_FUNCTION(Sin, sin, std::sin)
UNARY_FUNCTION(Cos, cos, std::cos)
UNARY_FUNCTION(Tan, tan, std::tan)
UNARY_FUNCTION(Asin, asin, std::asin)
UNARY_FUNCTION(Acos, acos, std::acos)
UNARY_FUNCTION(Atan, atan, std::atan)
UNARY_FUNCTION(Gamma, gamma, gsl_sf_gamma)

#undef UNARY_FUNCTION

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

EPtr UnaryOp::simplify() const {
    EPtr as = a->simplify();
    Constant* ac = dynamic_cast<Constant*>(as.get());
    if (ac != NULL) {
        return Constant::create(evaluate(ac->c));
    }
    return EPtr(construct(std::move(as)));
}

EPtr Neg::simplify() const {
    EPtr as = a->simplify();
    Constant* ac = dynamic_cast<Constant*>(as.get());
    if (ac != NULL) {
        return Constant::create(evaluate(ac->c));
    }
    Neg* an = dynamic_cast<Neg*>(as.get());
    if (an != NULL) {
        return an->a->ecopy();
    }
    return Neg::create(std::move(as));
}

EPtr Add::simplify() const {
    EPtr as = a->simplify();
    EPtr bs = b->simplify();
    Constant* ac = dynamic_cast<Constant*>(as.get());
    Constant* bc = dynamic_cast<Constant*>(bs.get());
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    if (ac && bc) {
        return Constant::create(acc + bcc);
    }
    // 0 + b = b
    if (ac && acc == 0) {
        return std::move(bs);
    }
    // a + 0 = a
    if (bc && bcc == 0) {
        return std::move(as);
    }
    Neg* an = dynamic_cast<Neg*>(as.get());
    Neg* bn = dynamic_cast<Neg*>(bs.get());
    // (-a) + (-b) = -(a+b)
    if (an && bn) {
        return Neg::create(Add::create(an->a->ecopy(), bn->a->ecopy()));
    }
    // (-a) + b = b-a
    if (an) {
        return Sub::create(std::move(bs), an->a->ecopy());
    }
    // a + (-b) = a-b
    if (bn) {
        return Sub::create(std::move(as), bn->a->ecopy());
    }
    return Add::create(std::move(as), std::move(bs));
}

EPtr Sub::simplify() const {
    std::unique_ptr<Add> _r(new Add(EPtr(a->copy()), EPtr(new Neg(EPtr(b->copy())))));
    return _r->simplify();
}

EPtr Mul::simplify() const {
    EPtr as(a->simplify());
    EPtr bs(b->simplify());
    Constant* ac = dynamic_cast<Constant*>(as.get());
    Constant* bc = dynamic_cast<Constant*>(bs.get());
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    // constant folding
    if (ac && bc) {
        return Constant::create(acc * bcc);
    }
    // 0 * b = 0
    // a * 0 = 0
    if ((ac && acc == 0) || (bc && bcc == 0)) {
        return EPtr(new Constant(0));
    }
    // 1 * b = b
    if (ac && acc == 1) {
        return std::move(bs);
    }
    // a * 1 = a
    if (bc && bcc == 1) {
        return std::move(as);
    }
    Neg* an = dynamic_cast<Neg*>(as.get());
    Neg* bn = dynamic_cast<Neg*>(bs.get());
    // (-a) * (-b) = a * b
    if (an && bn) {
        std::unique_ptr<Mul> _r(new Mul(EPtr(an->a->copy()), EPtr(bn->a->copy())));
        return _r->simplify();
    }
    // (-a) * b = -(a * b)
    if (an) {
        std::unique_ptr<Neg> _r(new Neg(EPtr(new Mul(EPtr(an->a->copy()), std::move(bs)))));
        return _r->simplify();
    }
    // a * (-b) = -(a * b)
    if (an) {
        std::unique_ptr<Neg> _r(new Neg(EPtr(new Mul(std::move(as), EPtr(bn->a->copy())))));
        return _r->simplify();
    }
    Div* ad = dynamic_cast<Div*>(as.get());
    Div* bd = dynamic_cast<Div*>(bs.get());
    // (aa/ab) * (ba/bb) = (aa * ba) / (ab * bb)
    if (ad && bd) {
        std::unique_ptr<Div> _r(new Div(EPtr(new Mul(EPtr(ad->a->copy()), EPtr(bd->a->copy()))), EPtr(new Mul(EPtr(ad->b->copy()), EPtr(bd->b->copy())))));
        return _r->simplify();
    }
    // (aa/ab) * b = (aa * b) / (ab)
    if (ad) {
        std::unique_ptr<Div> _r(new Div(EPtr(new Mul(EPtr(ad->a->copy()), std::move(bs))), EPtr(ad->b->copy())));
        return _r->simplify();
    }
    // a * (ba/bb) = (a * ba) / (bb)
    if (bd) {
        std::unique_ptr<Div> _r(new Div(EPtr(new Mul(std::move(as), EPtr(bd->a->copy()))), EPtr(bd->b->copy())));
        return _r->simplify();
    }  
    return EPtr(new Mul(std::move(as), std::move(bs)));
}

EPtr Div::simplify() const {
    EPtr as(a->simplify());
    EPtr bs(b->simplify());
    Constant* ac = dynamic_cast<Constant*>(as.get());
    Constant* bc = dynamic_cast<Constant*>(bs.get());
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    // constant folding
    if (ac && bc) {
        return EPtr(new Constant(acc / bcc));
    }
    // 0 / b = 0
    if (ac && acc == 0) {
        return EPtr(new Constant(0));
    }
    // a / 0 = NaN
    if (bc && bcc == 0) {
        return EPtr(new Constant(std::numeric_limits<Number>::quiet_NaN()));
    }
    // a / 1 = a
    if (bc && bcc == 1) {
        return std::move(as);
    }
    Neg* an = dynamic_cast<Neg*>(as.get());
    Neg* bn = dynamic_cast<Neg*>(bs.get());
    // (-a) / (-b) = a / b
    if (an && bn) {
        std::unique_ptr<Div> _r(new Div(EPtr(an->a->copy()), EPtr(bn->a->copy())));
        return _r->simplify();
    }
    // (-a) / b = -(a * b)
    if (an) {
        std::unique_ptr<Neg> _r(new Neg(EPtr(new Div(EPtr(an->a->copy()), std::move(bs)))));
        return _r->simplify();
    }
    // a * (-b) = -(a * b)
    if (an) {
        std::unique_ptr<Neg> _r(new Neg(EPtr(new Div(std::move(as), EPtr(bn->a->copy())))));
        return _r->simplify();
    }
    Div* ad = dynamic_cast<Div*>(as.get());
    Div* bd = dynamic_cast<Div*>(bs.get());
    // (aa/ab) / (ba/bb) = (aa * bb) / (ab * ba)
    if (ad && bd) {
        std::unique_ptr<Div> _r(new Div(Mul::create(ad->a->ecopy(), bd->b->ecopy()), Mul::create(ad->b->ecopy(), bd->a->ecopy())));
        return _r->simplify();
    }
    // (aa/ab) / b = aa / (ab * b)
    if (ad) {
        std::unique_ptr<Div> _r(new Div(EPtr(ad->a->copy()), EPtr(new Mul(EPtr(ad->b->copy()), std::move(bs)))));
        return _r->simplify();
    }
    // a / (ba/bb) = (a * bb) / (ba)
    if (bd) {
        std::unique_ptr<Div> _r(new Div(EPtr(new Mul(std::move(as), EPtr(bd->b->copy()))), EPtr(bd->a->copy())));
        return _r->simplify();
    }  
    return EPtr(new Div(std::move(as), std::move(bs)));
}

EPtr Pow::simplify() const {
    EPtr as(a->simplify());
    EPtr bs(b->simplify());
    Constant* ac = dynamic_cast<Constant*>(as.get());
    Constant* bc = dynamic_cast<Constant*>(bs.get());
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    // constant folding
    if (ac && bc) {
        return EPtr(new Constant(std::pow(acc, bcc)));
    }
    // a ^ 0 = 1; (0^0 = 1)
    // 1 ^ b = 1
    if ((bc && bcc == 0) || (ac && acc == 1)) {
        return EPtr(new Constant(1));
    }
/*    // 0 ^ b = 0
    if (ac && acc == 0) {
        return EPtr(new Constant(0));
    }*/
    // a ^ 1 = a
    if (bc && bcc == 1) {
        return std::move(as);
    }
    return EPtr(new Pow(std::move(as), std::move(bs)));
}

EPtr PowInt::simplify() const {
    EPtr as(a->simplify());
    Constant* ac = dynamic_cast<Constant*>(as.get());
    Number acc = ac ? ac->c : 0;
    // constant folding
    if (ac) {
        return EPtr(new Constant(powi(acc, b)));
    }
    // a ^ 0 = 1; (0^0 = 1)
    // 1 ^ b = 1
    if ((b == 0) || (ac && acc == 1)) {
        return EPtr(new Constant(1));
    }
/*    // 0 ^ b = 0
    if (ac && acc == 0) {
        return EPtr(new Constant(0));
    }*/
    // a ^ 1 = a
    if (b == 1) {
        return std::move(as);
    }
    return EPtr(new PowInt(std::move(as), b));
}

// d(-u)/dx = -du/dx
EPtr Neg::derivative(const Variable& var) const {
    return Neg::create(a->derivative(var));
}

// de^u/dx = de^u/du * du/dx = e^u * du/dx
EPtr Exp::derivative(const Variable& var) const {
    return Mul::create(ecopy(), a->derivative(var));
}
// dln(u)/dx = dln(u)/du * du/dx = 1/u * du/dx
EPtr Ln::derivative(const Variable& var) const {
    return Mul::create(Div::create(Constant::create(1), a->ecopy()), a->derivative(var));
}
// dsqrt(u)/dx = dsqrt(u)/du * du/dx = (1/2*sqrt(u)) * du/dx = (du/dx) / (2 * sqrt(u))
EPtr Sqrt::derivative(const Variable& var) const {
    return Div::create(a->derivative(var), Mul::create(Constant::create(2), Sqrt::create(a->ecopy())));
}
// dsin(u)/dx = dsin(u)/du * du/dx = cos(u) * du/dx
EPtr Sin::derivative(const Variable& var) const {
    return Mul::create(Cos::create(a->ecopy()), a->derivative(var));
}
// dcos(u)/dx = dcos(u)/du * du/dx = -sin(u) * du/dx
EPtr Cos::derivative(const Variable& var) const {
    return Neg::create(Mul::create(Sin::create(a->ecopy()), a->derivative(var)));
}
// dtan(u)/dx = dtan(u)/du * du/dx = (1 + tan^2(u)) * du/dx
EPtr Tan::derivative(const Variable& var) const {
    return Mul::create(Add::create(Constant::create(1), PowInt::create(Tan::create(a->ecopy()), 2)), a->derivative(var));
}
// d asin(u)/dx = d asin(u) / du * du/dx = (du/dx)/sqrt(1 - u^2)
EPtr Asin::derivative(const Variable& var) const {
    return Div::create(a->derivative(var), Sqrt::create(Sub::create(Constant::create(1), PowInt::create(a->ecopy(), 2))));
}
// d acos(u)/dx = d acos(u) / du * du/dx = -(du/dx)/sqrt(1 - u^2)
EPtr Acos::derivative(const Variable& var) const {
    return Neg::create(Div::create(a->derivative(var), Sqrt::create(Sub::create(Constant::create(1), PowInt::create(a->ecopy(), 2)))));
}
// d atan(u)/dx = d atan(u) / du * du/dx = (du/dx)/(1 + u^2)
EPtr Atan::derivative(const Variable& var) const {
    return Div::create(a->derivative(var), Add::create(Constant::create(1), PowInt::create(a->ecopy(), 2)));
}

// d(u+v)/dx = du/dx + dv/dx
EPtr Add::derivative(const Variable& var) const {
    return Add::create(a->derivative(var), b->derivative(var));
}
// d(u-v)/dx = du/dx - dv/dx
EPtr Sub::derivative(const Variable& var) const {
    return Sub::create(a->derivative(var), b->derivative(var));
}
// d(u*v)/dx = udv/dx + vdu/dx
EPtr Mul::derivative(const Variable& var) const {
    return Add::create(Mul::create(a->ecopy(), b->derivative(var)), Mul::create(b->ecopy(), a->derivative(var)));
}
// d(u/v)/dx = (vdu/dx - udv/dx) / (u^2)
EPtr Div::derivative(const Variable& var) const {
    return Div::create(
                Sub::create(Mul::create(b->ecopy(), a->derivative(var)), Mul::create(a->ecopy(), b->derivative(var))),
                PowInt::create(b->ecopy(), 2)
            );
}
// d(u^v)/dx = d/dx e^(vlogu) = d(e^(vlogu))/d(vlogu) d(vlogu)/dx
// = e^vlogu * (logu dv/dx + v dlogu/dx) = u^v * (logu dv/dx + v dlogu/du du/dx)
// = u^v * (logu dv/dx + v/u du/dx)
EPtr Pow::derivative(const Variable& var) const {
    return Mul::create(
                ecopy(),
                Add::create(
                    Mul::create(Ln::create(a->ecopy()), b->derivative(var)),
                    Mul::create(Div::create(b->ecopy(), a->ecopy()), a->derivative(var))
                )
            );
}
// d(u^a)/dx = d(u^b)/du * du/dx = bu^(b-1) * du/dx
EPtr PowInt::derivative(const Variable& var) const {
    return Mul::create(
                Mul::create(Constant::create(b), PowInt::create(a->ecopy(), b-1)),
                a->derivative(var)
            );
}

// dgamma(u)/dx = dgamma(u)/du * du/dx = gamma(u)*psi_0(u) * du/dx
EPtr Gamma::derivative(const Variable& var) const {
    return Mul::create(Mul::create(Gamma::create(a->ecopy()), PolyGamma::create(a->ecopy(), 0)), a->derivative(var));
}

// dpsi_b(u)/dx = dpsi_b(u)/du * du/dx = psi_(b+1)(u) * du/dx
EPtr PolyGamma::derivative(const Variable& var) const {
    return Mul::create(PolyGamma::create(a->ecopy(), b+1), a->derivative(var));
}
