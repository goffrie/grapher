#include "Expression.h"

#include <mmintrin.h>
#include <xmmintrin.h>

using namespace AsmJit;

static_assert(sizeof(sysint_t) == sizeof(void*), "sysint_t has the wrong size");

EvalFunc Expression::evaluator() const {
    Compiler c;
    c.setLogger(new FileLogger(stdout));
    c.newFunction(CALL_CONV_DEFAULT, FunctionBuilder1<float, int>()); // HACK: returning a "float" but really a v4sf
    c.ret(evaluate(c, c.argGP(0)));
    c.endFunction();
    return (EvalFunc)c.make();
}

XMMVar Variable::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = c.newXMM();
    switch (id->type) {
        case Id::Constant:
            c.movss(ret, dword_ptr_abs((void*) id->p));
            c.shufps(ret, ret, imm(0x00)); // replicate least significant data element
            return ret;
        case Id::Vector:
            c.movaps(ret, dqword_ptr_abs((void*) id->p, i, 4));
            return ret;
        default:
            throw this;
    }
}

XMMVar Constant::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = c.newXMM();
    c.movss(ret, dword_ptr_abs((void*) &c));
    c.shufps(ret, ret, imm(0x00)); // replicate least significant data element
    return ret;
}

XMMVar Neg::evaluate(Compiler& c, const GPVar& i) const {
    static const __v4si neg = {1 << 31, 1 << 31, 1 << 31, 1 << 31};
    XMMVar ret = a->evaluate(c, i);
    c.xorps(ret, dqword_ptr_abs((void*) &neg));
    return ret;
}

XMMVar Exp::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Ln::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Sqrt::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    c.sqrtps(ret, ret);
    return ret;
}

XMMVar Sin::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Cos::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Tan::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Asin::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Acos::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Atan::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Gamma::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Add::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i),
           xb = b->evaluate(c, i);
    c.addps(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Sub::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i),
           xb = b->evaluate(c, i);
    c.subps(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Mul::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i),
           xb = b->evaluate(c, i);
    c.mulps(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Div::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i),
           xb = b->evaluate(c, i);
    c.divps(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar PowInt::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar _a = a->evaluate(c, i),
           ret = c.newXMM();
    int _b = b;
    static const __v4sf ones = {1.f, 1.f, 1.f, 1.f};
    c.movaps(ret, dqword_ptr_abs((void*) &ones));
    bool neg = _b < 0;
    if (neg) _b = -_b;
    while (_b > 0) {
        if (_b & 1) c.mulps(ret, _a);
        _b >>= 1;
        if (_b > 0) c.mulps(_a, _a);
    }
    c.unuse(_a);
    if (neg) c.rcpps(ret, ret);
    return ret;
}

XMMVar Polynomial::evaluate(Compiler& c, const GPVar& i) const {
    if (!left.get()) return right->evaluate(c, i);
    XMMVar _l = left->evaluate(c, i);
    XMMVar _c = c.newXMM();
    switch (var.id->type) {
        case Variable::Id::Constant: {
            c.movss(_c, dword_ptr_abs((void*) var.id->p));
            c.shufps(_c, _c, imm(0x00));
            break;
        }
        case Variable::Id::Vector: {
            c.movaps(_c, dqword_ptr_abs((void*) var.id->p, i, 4));
            break;
        }
        default: throw var;
    }
    c.mulps(_l, _c);
    c.unuse(_c);
    XMMVar _r = right->evaluate(c, i);
    c.addps(_l, _r);
    c.unuse(_r);
    return _l;
}

XMMVar Pow::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar PolyGamma::evaluate(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluate(c, i);
    throw "Not implemented";
    return ret;
}
