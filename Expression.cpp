#include "Expression.h"

#include <cstring>

#include "util.h"
#include "global.h"
// #include "sse_mathfun.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>

//#define DEBUG_SIMPLIFY std::cerr << "simplifying at " << __LINE__ << " in " << __PRETTY_FUNCTION__ << ": " << toString() << std::endl;
#define DEBUG_SIMPLIFY

std::string tostr(float f) {
    std::ostringstream o;
    o << std::showpoint << std::setprecision(1) << f;
    return o.str();
}
std::string tostr(int i) {
    std::ostringstream o;
    o << i;
    return o.str();
}

std::string Constant::toString(int prec) const {
    return wrap(prec, (c <= -0) ? Precedence::Mul : Precedence::Const, tostr(c));
}

std::string PowInt::toString(int prec) const {
    return wrap(prec, Precedence::Pow, a->toString(Precedence::Pow+1) + " ^ " + tostr(b));
}

std::string PolyGamma::toString(int prec) const {
    return wrap(prec, Precedence::Func, std::string("psi_") + tostr(b) + "(" + a->toString(-1) + ")");
}



Vector Constant::evaluateVector(uz size) const {
    VectorR r = VECTOR_ALLOC(size);
    const Vc::float_v vc = c;
    VECTOR_LOOP V(r) = vc;
    return r;
}

Number Variable::evaluate() const {
    switch (id->type) {
        case Id::Constant:
        case Id::Vector:
            return *id->p;
        default:
            throw this;
    }
}

Vector Variable::evaluateVector(uz size) const {
    switch (id->type) {
        case Id::Constant: {
            VectorR r = VECTOR_ALLOC(size);
            const Vc::float_v vc = *id->p;
            VECTOR_LOOP V(r) = vc;
            return r;
        }
        case Id::Vector: {
            VectorR r = VECTOR_ALLOC(size);
            std::memcpy(r, id->p, size*sizeof(Number));
            return r;
        }
        default:
            throw this;
    }
}


#define UNARY_EVALUATE(Name, sfunc) \
Number Name::evaluate(Number _a) const { return sfunc(_a); } \
Number Name::evaluate() const { return sfunc(a->evaluate()); }
#define UNARY_FUNCTION(Name, sfunc) \
UNARY_EVALUATE(Name, sfunc) \
Vector Name::evaluateVector(uz size) const { \
    VectorR _a = a->evaluateVector(size); \
    for (uz i = 0; i < size; ++i) _a[i] = sfunc(_a[i]); \
    return _a; \
}

UNARY_EVALUATE(Neg, -)
UNARY_EVALUATE(Exp, std::exp)
UNARY_EVALUATE(Ln, std::log)
UNARY_EVALUATE(Sqrt, std::sqrt)
UNARY_EVALUATE(Sin, std::sin)
UNARY_EVALUATE(Cos, std::cos)
UNARY_EVALUATE(Tan, std::tan)
UNARY_EVALUATE(Asin, std::asin)
UNARY_EVALUATE(Acos, std::acos)
UNARY_EVALUATE(Atan, std::atan)
UNARY_FUNCTION(Gamma, gsl_sf_gamma)

#undef UNARY_FUNCTION
#undef UNARY_EVALUATE

#define UNARY_VECTOR_EVALUATE(Name, vfunc) \
Vector Name::evaluateVector(uz size) const { \
    VectorR _a = a->evaluateVector(size); \
    VECTOR_LOOP V(_a) = vfunc(V(_a)); \
    return _a; \
}

Vc::float_v tan(Vc::float_v v) {
    return Vc::sin(v) / Vc::cos(v);
}
Vc::float_v acos(Vc::float_v v) {
    static const Vc::float_v halfpi(M_PI/2.);
    return halfpi - Vc::asin(v);
}

UNARY_VECTOR_EVALUATE(Neg, -)
UNARY_VECTOR_EVALUATE(Sqrt, Vc::sqrt)
UNARY_VECTOR_EVALUATE(Ln, Vc::log)
UNARY_VECTOR_EVALUATE(Exp, Vc::exp)
UNARY_VECTOR_EVALUATE(Sin, Vc::sin)
UNARY_VECTOR_EVALUATE(Cos, Vc::cos)
UNARY_VECTOR_EVALUATE(Tan, tan)
UNARY_VECTOR_EVALUATE(Asin, Vc::asin)
UNARY_VECTOR_EVALUATE(Acos, acos)
UNARY_VECTOR_EVALUATE(Atan, Vc::atan)

#undef UNARY_VECTOR_EVALUATE

#define BINARY_VECTOR_EVALUATE(Name, op) \
Vector Name::evaluateVector(uz size) const { \
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
            VECTOR_LOOP reinterpret_cast<Vc::float_m&>(V(_a)) = V(_a) op V(_b); \
            break;
    switch (type) {
        IOP(LT, <)
        IOP(GT, >)
        IOP(LTE, <=)
        IOP(GTE, >=)
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

Vector PowInt::evaluateVector(uz size) const {
    VectorR _a = a->evaluateVector(size);
    int _b = b;
    VectorR ret = VECTOR_ALLOC(size);
    const Vc::float_v ones(Vc::One);
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

Vector Pow::evaluateVector(uz size) const {
    VectorR _a = a->evaluateVector(size);
    VectorR _b = b->evaluateVector(size);
    //for (uz i = 0; i < size; ++i) _a[i] = std::pow(_a[i], _b[i]);
    VECTOR_LOOP V(_a) = Vc::exp(V(_b) * Vc::log(V(_a)));
    VECTOR_FREE(_b);
    return _a;
}

Number PolyGamma::evaluate() const {
    return gsl_sf_psi_n(b, a->evaluate());
}

Number PolyGamma::evaluate(Number _a) const {
    return gsl_sf_psi_n(b, _a);
}

Vector PolyGamma::evaluateVector(uz size) const {
    VectorR _a = a->evaluateVector(size);
    for (uz i = 0; i < size; ++i) _a[i] = gsl_sf_psi_n(b, _a[i]);
    return _a;
}

EPtr UnaryOp::simplify() const {
    DEBUG_SIMPLIFY
    EPtr as = a->simplify();
    Constant* ac = dynamic_cast<Constant*>(as.get());
    if (ac != NULL) {
        return Constant::create(evaluate(ac->c));
    }
    return EPtr(construct(std::move(as)));
}

EPtr Neg::simplify() const {
    DEBUG_SIMPLIFY
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
    DEBUG_SIMPLIFY
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
    DEBUG_SIMPLIFY
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
    DEBUG_SIMPLIFY
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
    if (bn) {
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
    
    Mul* am = dynamic_cast<Mul*>(as.get());
    Mul* bm = dynamic_cast<Mul*>(bs.get());
    if (am && bm) {
        Constant* amc = dynamic_cast<Constant*>(am->a.get());
        Constant* bmc = dynamic_cast<Constant*>(bm->a.get());
        // (a * x) * (b * y) = (a * b) * x * y
        if (amc && bmc) {
            return (Constant::create(amc->c * bmc->c) * (std::move(am->b) * std::move(bm->b)))->simplify();
        }
        if (bmc) { // move the constant to the left
            return (std::move(bm->a) * std::move(as) * std::move(bm->b))->simplify();
        }
    }
    if (bc && !ac) return std::move(bs) * std::move(as); // put the constant on the left
    return std::move(as) * std::move(bs);
}

EPtr Div::simplify() const {
    DEBUG_SIMPLIFY
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
    DEBUG_SIMPLIFY
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
        return PowInt::create(std::move(as), qRound(bcc))->simplify();
    }
    return Pow::create(std::move(as), std::move(bs));
}

EPtr PowInt::simplify() const {
    DEBUG_SIMPLIFY
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
    return (PowInt::create(Cos::create(a->ecopy()), -2)) * a->derivative(var);
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

bool Pow::polynomial(const Variable& var) const {
    return a->polynomial(var) && dynamic_cast<Constant*>(b.get()); // maybe too strict
}

bool Div::polynomial(const Variable& var) const {
    Constant* bc = dynamic_cast<Constant*>(b.get());
    Variable* bv = dynamic_cast<Variable*>(b.get());
    return a->polynomial(var) && (bc || (bv && *bv != var)); // perhaps try polynomial division
}

EPtr Neg::expand() const {
    EPtr _a = a->expand();
    Add* aa = dynamic_cast<Add*>(_a.get());
    if (aa) {
        return ((-std::move(aa->a)) + (-std::move(aa->b)))->expand();
    }
    return -std::move(_a);
}

EPtr PowInt::expand() const {
    if (b < 0) return Constant::create(1) / (PowInt::create(a->ecopy(), -b))->expand();
    if (b == 0) return Constant::create(1);
    EPtr _a = a->expand();
    if (b == 1) return _a;
    Add* aa = dynamic_cast<Add*>(_a.get());
    if (aa) {
        EPtr ret = PowInt::create(aa->a->ecopy(), b);
        for (int i = 1; i < b; ++i) {
            ret = std::move(ret) + Constant::create(binom(b, i)) * PowInt::create(aa->a->ecopy(), b-i) * PowInt::create(aa->b->ecopy(), i);
        }
        ret = std::move(ret) + PowInt::create(aa->b->ecopy(), b);
        return ret->expand();
    }
    Mul* am = dynamic_cast<Mul*>(_a.get());
    if (am) {
        return (PowInt::create(std::move(am->a), b) * PowInt::create(std::move(am->b), b))->expand();
    }
    Div* ad = dynamic_cast<Div*>(_a.get());
    if (ad) {
        return (PowInt::create(std::move(am->a), b) / PowInt::create(std::move(am->b), b))->expand();
    }
    Neg* an = dynamic_cast<Neg*>(_a.get());
    if (an) {
        if (b & 1) {
            return -(PowInt::create(std::move(an->a), b)->expand());
        } else {
            return PowInt::create(std::move(an->a), b)->expand();
        }
    }
    PowInt* api = dynamic_cast<PowInt*>(_a.get());
    if (api) {
        return PowInt::create(std::move(api->a), b*api->b)->expand();
    }
    Pow* ap = dynamic_cast<Pow*>(_a.get());
    if (ap) {
        return Pow::create(std::move(ap->a), Constant::create(b)*std::move(ap->b))->expand();
    }
    return PowInt::create(std::move(_a), b);
}

EPtr Add::expand() const {
    return a->expand() + b->expand();
}

EPtr Sub::expand() const {
    return a->expand() + (-b->ecopy())->expand();
}

EPtr Div::expand() const {
    EPtr _a = a->expand();
    EPtr _b = b->expand();
    Constant* bc = dynamic_cast<Constant*>(_b.get());
    Variable* bv = dynamic_cast<Variable*>(_b.get());
    Add* aa = dynamic_cast<Add*>(_a.get());
    if (((bc || bv) && aa)) {
        EPtr _bcopy = _b->ecopy();
        return ((std::move(aa->a) / std::move(_b)) + (std::move(aa->b) / std::move(_bcopy)))->expand();
    } else {
        return std::move(_a) / std::move(_b);
    }
}

EPtr Mul::expand() const {
    EPtr _a = a->expand();
    EPtr _b = b->expand();
    Add* a_add = dynamic_cast<Add*>(_a.get());
    Add* b_add = dynamic_cast<Add*>(_b.get());
    if (a_add && b_add) {
        EPtr& aa = a_add->a;
        EPtr& ab = a_add->b;
        EPtr& ba = b_add->a;
        EPtr& bb = b_add->b;
        EPtr aa2 = aa->ecopy();
        EPtr ab2 = ab->ecopy();
        EPtr ba2 = ba->ecopy();
        EPtr bb2 = bb->ecopy();
        return (std::move(aa)*std::move(ba)+std::move(aa2)*std::move(bb)+std::move(ab)*std::move(ba2)+std::move(ab2)*std::move(bb2))->expand();
    }
    if (a_add) {
        EPtr& aa = a_add->a;
        EPtr& ab = a_add->b;
        EPtr b2 = _b->ecopy();
        return (std::move(aa)*std::move(_b)+std::move(ab)*std::move(b2))->expand();
    }
    if (b_add) {
        EPtr& ba = b_add->a;
        EPtr& bb = b_add->b;
        EPtr a2 = _a->ecopy();
        return (std::move(ba)*std::move(_a)+std::move(bb)*std::move(a2))->expand();
    }
    return std::move(_a) * std::move(_b);
}

std::unique_ptr<Polynomial> Expression::facsum(const Variable& var) const {
    return Polynomial::create(var, simplify());
}

std::unique_ptr<Polynomial> Add::facsum(const Variable& var) const {
    std::cerr << "ENTER Add::facsum(" << this << "): " << toString() << std::endl;
    std::unique_ptr<Polynomial> _a = a->facsum(var);
    std::unique_ptr<Polynomial> _b = b->facsum(var);
    std::cerr << "Add::facsum: (" << _a->toString() << ") + (" << _b->toString() << ")" << std::endl;
    if (_b->left.get()) {
        if (_a->left.get()) {
            _a->left = (std::move(_a->left) + std::move(_b->left))->facsum(var);
        } else {
            _a->left = std::move(_b->left);
        }
    }
    _a->right = (std::move(_a->right) + std::move(_b->right))->simplify();
    std::cerr << "LEAVE Add::facsum(" << this << "): " << toString() << " -> " << _a->toString() << std::endl;
    return std::move(_a);
}
std::unique_ptr<Polynomial> Mul::facsum(const Variable& var) const {
    std::cerr << "ENTER Mul::facsum(" << this << "): " << toString() << std::endl;
    std::unique_ptr<Polynomial> _a = a->facsum(var);
    std::unique_ptr<Polynomial> _b = b->facsum(var);
    std::cerr << "Mul::facsum: (" << _a->toString() << ") * (" << _b->toString() << ")" << std::endl;
    std::unique_ptr<Polynomial> ret;
    if (_a->left.get() && _b->left.get()) {
        EPtr aa = std::move(_a->left);
        EPtr ab = std::move(_a->right);
        EPtr ba = std::move(_b->left);
        EPtr bb = std::move(_b->right);
        // (aa*x+ab)(ba*x+bb) = (aa*ba)*x^2+(aa*bb+ab*ba)*x+(ab*bb)
        std::cerr << "Mul::facsum(" << this << "): mid: ";
        std::unique_ptr<Polynomial> mid = (aa->ecopy()*bb->ecopy()+ab->ecopy()*ba->ecopy())->facsum(var); // (aa*bb+ab*ba) = mid->left * x + mid->right
        std::cerr << "Mul::facsum(" << this << "): left: ";
        std::unique_ptr<Polynomial> left = (std::move(aa)*std::move(ba))->facsum(var); // aa*ba = left->left * x + left->right
        if (mid->left.get()) {
            left = (std::move(left) + std::move(mid->left))->facsum(var);
            // nleft = left->left * x + (left->right + mid->left)
        }
        EPtr right = (std::move(ab)*std::move(bb))->simplify(); // (ab * bb) = right
        // (aa*ba)x^2 + (aa*bb + ab*ba)*x + (ab*bb)
        // = (left->left * x + left->right) * x^2 + (mid->left * x + mid->right) * x + right
        // = left->left * x^3 + left->right * x^2 + mid->left * x^2 + mid->right * x + right
        // = left->left * x^3 + (left->right + mid->left) * x^2 + mid->right * x + right
        // = ((left->left * x + (left->right + mid->left)) * x + mid->right) * x + right
        // = (nleft * x + mid->right) * x + right
        ret = Polynomial::create(Polynomial::create(std::move(left), var, std::move(mid->right)), var, std::move(right));
    } else if (_a->left.get()) {
        // (aa*x+ab) * b = (aa * b) * x + (ab * b)
        EPtr _b2 = _b->ecopy();
        ret = Polynomial::create((std::move(_a->left) * std::move(_b))->facsum(var), var, (std::move(_a->right) * std::move(_b2))->simplify());
    } else if (_b->left.get()) {
        // (ba*x+bb) * a = (ba * a) * x + (bb * a)
        EPtr _a2 = _a->ecopy();
        ret = Polynomial::create((std::move(_b->left) * std::move(_a))->facsum(var), var, (std::move(_b->right) * std::move(_a2))->simplify());
    } else {
        ret = Polynomial::create(var, (std::move(_a) * std::move(_b))->simplify());
    }
    std::cerr << "LEAVE Mul::facsum(" << this << "): " << toString() << " -> " << ret->toString() << std::endl;
    return std::move(ret);
}
std::unique_ptr<Polynomial> Div::facsum(const Variable& var) const {
    std::unique_ptr<Polynomial> _a = a->facsum(var);
    std::unique_ptr<Polynomial> _b = b->facsum(var);
    if (_a->left.get()) {
        if (_b->left.get()) {
            // dang
            return Polynomial::create(var, std::move(_a) / std::move(_b));
        } else {
            // (aa * x + ab) / bb = (aa / bb) * x + (ab / bb)
            EPtr _b2 = _b->right->ecopy();
            return Polynomial::create((std::move(_a->left) / std::move(_b->right))->facsum(var), var, (std::move(_a->right) / std::move(_b2))->simplify());
        }
    } else {
        if (_b->left.get()) {
            // dangdang
            return Polynomial::create(var, (std::move(_a->right) / std::move(_b))->simplify());
        } else {
            return Polynomial::create(var, (std::move(_a->right) / std::move(_b->right))->simplify());
        }
    }
}

std::unique_ptr<Polynomial> Neg::facsum(const Variable& var) const {
    std::unique_ptr<Polynomial> _a = a->facsum(var);
    if (_a->left.get()) {
        return Polynomial::create((-std::move(_a->left))->facsum(var), var, (-std::move(_a->right))->simplify());
    } else {
        return Polynomial::create(var, (-std::move(_a->right))->simplify());
    }
}

std::unique_ptr<Polynomial> Variable::facsum(const Variable& var) const {
    if (*this == var) {
        return Polynomial::create(Polynomial::create(var, Constant::create(1)), var, Constant::create(0));
    } else {
        return Polynomial::create(var, ecopy());
    }
}

std::unique_ptr<Polynomial> Polynomial::facsum(const Variable& _var) const {
    if (var == _var) {
        return std::unique_ptr<Polynomial>(copy());
    } else {
        return Polynomial::create(_var, ecopy());
    }
}

std::unique_ptr<Polynomial> PowInt::facsum(const Variable& var) const {
    if (b == 0) return Polynomial::create(var, Constant::create(1));
    if (b == 1) return a->facsum(var);
    Variable* _a = dynamic_cast<Variable*>(a.get());
    if (b < 0 || !_a || *_a != var) return Polynomial::create(var, ecopy());
    std::unique_ptr<Polynomial> r = _a->facsum(var);
    for (int i = 1; i < b; ++i) r = Polynomial::create(std::move(r), var, Constant::create(0));
    return std::move(r);
}

Polynomial* Polynomial::copy() const {
    if (left.get())
        return new Polynomial(std::unique_ptr<Polynomial>(left->copy()), var, right->ecopy());
    else
        return new Polynomial(var, right->ecopy());
}

EPtr Polynomial::derivative(const Variable& res) const {
    if (!left.get()) {
        if (res == var)
            return Polynomial::create(var, right->derivative(var));
        else
            return right->derivative(res);
    }
    if (res != var) {
        // d/dx (a * u + b) = da/dx * u + a * du/dx + db/dx
        //     [du/dx == 0]
        // = da/dx *u + db/dx
        return Polynomial::create(left->derivative(res), var, right->derivative(res));
    } else {
        // d/dx (ax+b) = da/dx*x+dx/dx*a = da/dx*x+a
        return Polynomial::create(left->derivative(var), var, left->ecopy());
    }
}

Number Polynomial::evaluate() const {
    if (!left.get()) return right->evaluate();
    return left->evaluate() * var.evaluate() + right->evaluate();
}

Vector Polynomial::evaluateVector(std::size_t size) const {
    if (!left.get()) return right->evaluateVector(size);
    VectorR _l = left->evaluateVector(size);
    switch (var.id->type) {
        case Variable::Id::Constant: {
            Vc::float_v c = *var.id->p;
            VECTOR_LOOP V(_l) *= c;
            break;
        }
        case Variable::Id::Vector: {
            VectorR c = var.id->p;
            VECTOR_LOOP V(_l) *= V(c);
            break;
        }
        default: throw var;
    }
    VectorR _r = right->evaluateVector(size);
    VECTOR_LOOP V(_l) += V(_r);
    VECTOR_FREE(_r);
    return _l;
}

bool Polynomial::polynomial(const Variable& v) const {
    if (left.get() && !left->polynomial(v)) return false;
    return right->polynomial(v);
}

EPtr Polynomial::substitute(const Expression::Subst& s) const {
    if (!left.get()) return Polynomial::create(var, right->substitute(s));
    EPtr _v = var.substitute(s);
    Variable* _vv = dynamic_cast<Variable*>(_v.get());
    if (!_vv || *_vv != var) {
        return left->substitute(s) * std::move(_v) + right->substitute(s);
    } else {
        return Polynomial::create(left->substitute(s), var, right->substitute(s));
    }
}

void Polynomial::variables(std::set<Variable>& out) const {
    if (left.get()) {
        left->variables(out);
        var.variables(out);
    }
    right->variables(out);
}

int Polynomial::degree() const {
    return left.get() ? (left->degree() + 1) : 0;
}

EPtr Polynomial::create(EPtr l, const Variable& v, EPtr r) {
    Polynomial* lp = dynamic_cast<Polynomial*>(l.get());
    if (lp) {
        l.release();
        return Polynomial::create(std::unique_ptr<Polynomial>(lp), v, std::move(r));
    } else {
        return std::move(l) * v.ecopy() + std::move(r);
    }
}

std::unique_ptr<Polynomial> Polynomial::create(std::unique_ptr<Polynomial> l, const Variable& v, EPtr r) {
    return std::unique_ptr<Polynomial>(new Polynomial(std::move(l), v, std::move(r)));
}

std::unique_ptr<Polynomial> Polynomial::create(const Variable& v, EPtr r) {
    return std::unique_ptr<Polynomial>(new Polynomial(v, std::move(r)));
}

EPtr Polynomial::simplify() const {
    DEBUG_SIMPLIFY
    if (!left.get()) return right->simplify();
    return Expression::simplify();
}

std::string Polynomial::toString(int prec) const {
    if (left.get()) {
        return wrap(prec, Precedence::Add, left->toString(Precedence::Mul) + " * " + var.toString(Precedence::Mul) + " + " + right->toString(Precedence::Add));
    } else {
        return right->toString(prec);
    }
}
