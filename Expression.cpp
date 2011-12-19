#include "Expression.h"

#include <cstring>

#include "util.h"

#include <xmmintrin.h>
#define USE_SSE2
#include "sse_mathfun.h"

EPtr operator-(EPtr a) {
    return Neg::create(std::move(a));
}
EPtr operator+(EPtr a, EPtr b) {
    return Add::create(std::move(a), std::move(b));
}
EPtr operator-(EPtr a, EPtr b) {
    return Sub::create(std::move(a), std::move(b));
}
EPtr operator*(EPtr a, EPtr b) {
    return Mul::create(std::move(a), std::move(b));
}
EPtr operator/(EPtr a, EPtr b) {
    return Div::create(std::move(a), std::move(b));
}

#define VECTOR_LOOP for (std::size_t i = 0; i < size; i += SSE_VECTOR_SIZE)
#define V(a) (*reinterpret_cast<v4sf*>(a+i))

Vector Constant::evaluateVector(size_t size) const {
    VectorR r = VECTOR_ALLOC(size);
    const v4sf _c = {c, c, c, c};
    VECTOR_LOOP V(r) = _c;
    return r;
}

Vector External::evaluateVector(size_t size) const {
    VectorR r = VECTOR_ALLOC(size);
    std::memcpy(r, c, sizeof(Number)*size);
    return r;
}

#define UNARY_EVALUATE(Name, sfunc) \
Number Name::evaluate(Number _a) const { return sfunc(_a); } \
Number Name::evaluate() const { return sfunc(a->evaluate()); }
#define UNARY_FUNCTION(Name, sfunc) \
UNARY_EVALUATE(Name, sfunc) \
Vector Name::evaluateVector(std::size_t size) const { \
    VectorR _a = a->evaluateVector(size); \
    for (std::size_t i = 0; i < size; ++i) _a[i] = sfunc(_a[i]); \
    return _a; \
}

UNARY_EVALUATE(Neg, -)
UNARY_EVALUATE(Exp, std::exp)
UNARY_EVALUATE(Ln, std::log)
UNARY_EVALUATE(Sqrt, std::sqrt)
UNARY_EVALUATE(Sin, std::sin)
UNARY_EVALUATE(Cos, std::cos)
UNARY_EVALUATE(Tan, std::tan)
UNARY_FUNCTION(Asin, std::asin)
UNARY_FUNCTION(Acos, std::acos)
UNARY_FUNCTION(Atan, std::atan)
UNARY_FUNCTION(Gamma, gsl_sf_gamma)

#undef UNARY_FUNCTION
#undef UNARY_EVALUATE

#define UNARY_VECTOR_EVALUATE(Name, vfunc) \
Vector Name::evaluateVector(std::size_t size) const { \
    VectorR _a = a->evaluateVector(size); \
    VECTOR_LOOP V(_a) = vfunc(V(_a)); \
    return _a; \
}

UNARY_VECTOR_EVALUATE(Neg, -)
UNARY_VECTOR_EVALUATE(Sqrt, _mm_sqrt_ps)
UNARY_VECTOR_EVALUATE(Ln, log_ps)
UNARY_VECTOR_EVALUATE(Exp, exp_ps)
UNARY_VECTOR_EVALUATE(Sin, sin_ps)
UNARY_VECTOR_EVALUATE(Cos, cos_ps)
inline v4sf tan_ps(v4sf x) {
    v4sf sin, cos;
    sincos_ps(x, &sin, &cos);
    return sin / cos;
}
UNARY_VECTOR_EVALUATE(Tan, tan_ps)

#undef UNARY_VECTOR_EVALUATE

#define BINARY_VECTOR_EVALUATE(Name, op) \
Vector Name::evaluateVector(std::size_t size) const { \
    VectorR _a = a->evaluateVector(size); \
    VectorR _b = b->evaluateVector(size); \
    VECTOR_LOOP V(_a) op##= V(_b); \
    VECTOR_FREE(_b); \
    return _a; \
}

BINARY_VECTOR_EVALUATE(Add, +)
BINARY_VECTOR_EVALUATE(Sub, -)
BINARY_VECTOR_EVALUATE(Mul, *)
BINARY_VECTOR_EVALUATE(Div, /)

#undef BINARY_VECTOR_EVALUATE

Vector Inequality::evaluateVector(size_t size) const {
    VectorR _a = a->evaluateVector(size);
    VectorR _b = b->evaluateVector(size);
#define IOP(n, op) \
        case n: \
            VECTOR_LOOP V(_a) = _mm_cmp##op##_ps(V(_a), V(_b)); \
            break;
    switch (type) {
        IOP(LT, lt)
        IOP(GT, gt)
        IOP(LTE, le)
        IOP(GTE, ge)
    }
#undef IOP
    VECTOR_FREE(_b);
    return _a;
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
    VectorR ret = VECTOR_ALLOC(size);
    constexpr v4sf ones = {1.f, 1.f, 1.f, 1.f};
    VECTOR_LOOP V(ret) = ones;
    bool neg = _b < 0;
    if (neg) _b = -_b;
    while (_b > 0) {
        if (_b & 1) VECTOR_LOOP V(ret) *= V(_a);
        _b >>= 1;
        VECTOR_LOOP V(_a) *= V(_a);
    }
    if (neg) VECTOR_LOOP V(ret) = ones / V(ret);
    VECTOR_FREE(_a);
    return ret;
}

Vector Pow::evaluateVector(std::size_t size) const {
    VectorR _a = a->evaluateVector(size);
    VectorR _b = b->evaluateVector(size);
    for (std::size_t i = 0; i < size; ++i) _a[i] = std::pow(_a[i], _b[i]);
    VECTOR_FREE(_b);
    return _a;
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
        return std::move(an->a);
    }
    Sub* asb = dynamic_cast<Sub*>(as.get());
    if (asb != NULL) {
        return (std::move(asb->b) - std::move(asb->a))->simplify();
    }
    return -std::move(as);
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
        return (-(std::move(an->a) + std::move(bn->a)))->simplify();
    }
    // (-a) + b = b-a
    if (an) {
        return (std::move(bs) - std::move(an->a))->simplify();
    }
    // a + (-b) = a-b
    if (bn) {
        return (std::move(as) - std::move(bn->a))->simplify();
    }
    return std::move(as) + std::move(bs);
}

EPtr Sub::simplify() const {
    EPtr as = a->simplify();
    EPtr bs = b->simplify();
    Constant* ac = dynamic_cast<Constant*>(as.get());
    Constant* bc = dynamic_cast<Constant*>(bs.get());
    Number acc = ac ? ac->c : 0;
    Number bcc = bc ? bc->c : 0;
    if (ac && bc) {
        return Constant::create(acc - bcc);
    }
    // 0 - b = -b
    if (ac && acc == 0) {
        return (-std::move(bs))->simplify();
    }
    // a - 0 = a
    if (bc && bcc == 0) {
        return std::move(as);
    }
    Neg* an = dynamic_cast<Neg*>(as.get());
    Neg* bn = dynamic_cast<Neg*>(bs.get());
    // (-a) - (-b) = b-a
    if (an && bn) {
        return (std::move(bn->a) - std::move(an->a))->simplify();
    }
    // (-a) - b = -(a+b)
    if (an) {
        return (-(std::move(an->a) + std::move(bs)))->simplify();
    }
    // a - (-b) = a+b
    if (bn) {
        return (std::move(as) + std::move(bn->a))->simplify();
    }
    return std::move(as) - std::move(bs);
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
        return Constant::create(0);
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
        return (std::move(an->a) * std::move(bn->a))->simplify();
    }
    // (-a) * b = -(a * b)
    if (an) {
        return (-(std::move(an->a) * std::move(bs)))->simplify();
    }
    // a * (-b) = -(a * b)
    if (an) {
        return (-(std::move(as) * std::move(bn->a)))->simplify();
    }
    Div* ad = dynamic_cast<Div*>(as.get());
    Div* bd = dynamic_cast<Div*>(bs.get());
    // (aa/ab) * (ba/bb) = (aa * ba) / (ab * bb)
    if (ad && bd) {
        return ((std::move(ad->a) * std::move(bd->a)) / (std::move(ad->b) * std::move(bd->b)))->simplify();
    }
    // (aa/ab) * b = (aa * b) / (ab)
    if (ad) {
        return ((std::move(ad->a) * std::move(bs)) / std::move(ad->b))->simplify();
    }
    // a * (ba/bb) = (a * ba) / (bb)
    if (bd) {
        return ((std::move(as) * std::move(bd->a)) / std::move(bd->b))->simplify();
    }  
    return std::move(as) * std::move(bs);
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
        return Constant::create(acc / bcc);
    }
    // 0 / b = 0
    if (ac && acc == 0) {
        return Constant::create(0);
    }
    // a / 0 = +inf
    if (bc && bcc == 0) {
        return Constant::create(std::numeric_limits<Number>::infinity());
    }
    // a / 1 = a
    if (bc && bcc == 1) {
        return std::move(as);
    }
    Neg* an = dynamic_cast<Neg*>(as.get());
    Neg* bn = dynamic_cast<Neg*>(bs.get());
    // (-a) / (-b) = a / b
    if (an && bn) {
        return (std::move(an->a) / std::move(bn->a))->simplify();
    }
    // (-a) / b = -(a * b)
    if (an) {
        return (-(std::move(an->a) / std::move(bs)))->simplify();
    }
    // a * (-b) = -(a * b)
    if (an) {
        return (-(std::move(as) / std::move(bn->a)))->simplify();
    }
    Div* ad = dynamic_cast<Div*>(as.get());
    Div* bd = dynamic_cast<Div*>(bs.get());
    // (aa/ab) / (ba/bb) = (aa * bb) / (ab * ba)
    if (ad && bd) {
        return ((std::move(ad->a) * std::move(bd->b)) / (std::move(ad->b) * std::move(bd->a)))->simplify();
    }
    // (aa/ab) / b = aa / (ab * b)
    if (ad) {
        return (std::move(ad->a) / (std::move(ad->b) * std::move(bs)))->simplify();
    }
    // a / (ba/bb) = (a * bb) / (ba)
    if (bd) {
        return ((std::move(as) * std::move(bd->b)) / std::move(bd->a))->simplify();
    }  
    return std::move(as) / std::move(bs);
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
        return Constant::create(std::pow(acc, bcc));
    }
    // a ^ 0 = 1; (0^0 = 1)
    // 1 ^ b = 1
    if ((bc && bcc == 0) || (ac && acc == 1)) {
        return Constant::create(1);
    }
/*    // 0 ^ b = 0
    if (ac && acc == 0) {
        return Constant::create(0);
    }*/
    // a ^ 1 = a
    if (bc && bcc == 1) {
        return std::move(as);
    }
    if (bc && isIntegral(bcc)) {
        return PowInt::create(std::move(as), rnd(bcc))->simplify();
    }
    return Pow::create(std::move(as), std::move(bs));
}

EPtr PowInt::simplify() const {
    EPtr as(a->simplify());
    Constant* ac = dynamic_cast<Constant*>(as.get());
    Number acc = ac ? ac->c : 0;
    // constant folding
    if (ac) {
        return Constant::create(powi(acc, b));
    }
    // a ^ 0 = 1; (0^0 = 1)
    // 1 ^ b = 1
    if ((b == 0) || (ac && acc == 1)) {
        return Constant::create(1);
    }
/*    // 0 ^ b = 0
    if (ac && acc == 0) {
        return EPtr(new Constant(0));
    }*/
    // a ^ 1 = a
    if (b == 1) {
        return std::move(as);
    }
    return PowInt::create(std::move(as), b);
}

std::unique_ptr<Inequality> Inequality::simplify() const {
    return Inequality::create(Sub::create(a->simplify(), b->simplify())->simplify(), Constant::create(0), type);
}

// d(-u)/dx = -du/dx
EPtr Neg::derivative(const Variable& var) const {
    return -a->derivative(var);
}

// de^u/dx = de^u/du * du/dx = e^u * du/dx
EPtr Exp::derivative(const Variable& var) const {
    return ecopy() * a->derivative(var);
}
// dln(u)/dx = dln(u)/du * du/dx = 1/u * du/dx
EPtr Ln::derivative(const Variable& var) const {
    return (Constant::create(1) / a->ecopy()) * a->derivative(var);
}
// dsqrt(u)/dx = dsqrt(u)/du * du/dx = (1/2*sqrt(u)) * du/dx = (du/dx) / (2 * sqrt(u))
EPtr Sqrt::derivative(const Variable& var) const {
    return a->derivative(var) / (Constant::create(2) * Sqrt::create(a->ecopy()));
}
// dsin(u)/dx = dsin(u)/du * du/dx = cos(u) * du/dx
EPtr Sin::derivative(const Variable& var) const {
    return Cos::create(a->ecopy()) * a->derivative(var);
}
// dcos(u)/dx = dcos(u)/du * du/dx = -sin(u) * du/dx
EPtr Cos::derivative(const Variable& var) const {
    return -(Sin::create(a->ecopy()) * a->derivative(var));
}
// dtan(u)/dx = dtan(u)/du * du/dx = (1 + tan^2(u)) * du/dx
EPtr Tan::derivative(const Variable& var) const {
    return (Constant::create(1) + PowInt::create(Tan::create(a->ecopy()), 2)) * a->derivative(var);
}
// d asin(u)/dx = d asin(u) / du * du/dx = (du/dx)/sqrt(1 - u^2)
EPtr Asin::derivative(const Variable& var) const {
    return a->derivative(var) / Sqrt::create(Constant::create(1) - PowInt::create(a->ecopy(), 2));
}
// d acos(u)/dx = d acos(u) / du * du/dx = -(du/dx)/sqrt(1 - u^2)
EPtr Acos::derivative(const Variable& var) const {
    return -(a->derivative(var) / Sqrt::create(Constant::create(1) - PowInt::create(a->ecopy(), 2)));
}
// d atan(u)/dx = d atan(u) / du * du/dx = (du/dx)/(1 + u^2)
EPtr Atan::derivative(const Variable& var) const {
    return a->derivative(var) / (Constant::create(1) + PowInt::create(a->ecopy(), 2));
}

// d(u+v)/dx = du/dx + dv/dx
EPtr Add::derivative(const Variable& var) const {
    return a->derivative(var) + b->derivative(var);
}
// d(u-v)/dx = du/dx - dv/dx
EPtr Sub::derivative(const Variable& var) const {
    return a->derivative(var) - b->derivative(var);
}
// d(u*v)/dx = udv/dx + vdu/dx
EPtr Mul::derivative(const Variable& var) const {
    return a->ecopy() * b->derivative(var) + b->ecopy() * a->derivative(var);
}
// d(u/v)/dx = (vdu/dx - udv/dx) / (u^2)
EPtr Div::derivative(const Variable& var) const {
    return (b->ecopy() * a->derivative(var) - a->ecopy() * b->derivative(var)) / PowInt::create(b->ecopy(), 2);
}
// d(u^v)/dx = d/dx e^(vlogu) = d(e^(vlogu))/d(vlogu) d(vlogu)/dx
// = e^vlogu * (logu dv/dx + v dlogu/dx) = u^v * (logu dv/dx + v dlogu/du du/dx)
// = u^v * (logu dv/dx + v/u du/dx)
EPtr Pow::derivative(const Variable& var) const {
    return ecopy() * (Ln::create(a->ecopy()) * b->derivative(var) + b->ecopy() / a->ecopy() * a->derivative(var));
}
// d(u^a)/dx = d(u^b)/du * du/dx = b * u^(b-1) * du/dx
EPtr PowInt::derivative(const Variable& var) const {
    return Constant::create(b) * PowInt::create(a->ecopy(), b-1) * a->derivative(var);
}

// dgamma(u)/dx = dgamma(u)/du * du/dx = gamma(u)*psi_0(u) * du/dx
EPtr Gamma::derivative(const Variable& var) const {
    return Gamma::create(a->ecopy()) * PolyGamma::create(a->ecopy(), 0) * a->derivative(var);
}

// dpsi_b(u)/dx = dpsi_b(u)/du * du/dx = psi_(b+1)(u) * du/dx
EPtr PolyGamma::derivative(const Variable& var) const {
    return PolyGamma::create(a->ecopy(), b+1) * a->derivative(var);
}
