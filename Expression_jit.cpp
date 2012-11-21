#include "Expression.h"

using namespace AsmJit;

static_assert(sizeof(sysint_t) == sizeof(void*), "sysint_t has the wrong size");

EvalFunc Expression::evaluator() const {
    Compiler c;
    c.setLogger(new FileLogger(stdout));
    c.newFunction(CALL_CONV_DEFAULT, FunctionBuilder0<double>());
    c.ret(evaluate(c));
    c.endFunction();
    return (EvalFunc)c.make();
}

XMMVar Variable::evaluate(Compiler& c) const {
    XMMVar ret = c.newXMM();
    switch (id->type) {
        case Id::Constant:
        case Id::Vector:
#ifdef ASMJIT_X64
// stupid x86-64 quirk: no absolute 64-bit addressing for most instructions
            if ((sysuint_t) id->p > (sysuint_t) 0xFFFFFFFF) {
                GPVar temp = c.newGP();
                c.mov(temp, imm((sysint_t) (void*) id->p));
                c.movss(ret, dword_ptr(temp));
                c.unuse(temp);
            } else
#endif
                c.movss(ret, dword_ptr_abs((void*) id->p));
            c.cvtss2sd(ret, ret); // convert float to double
            return ret;
        default:
            throw this;
    }
}

XMMVar Constant::evaluate(Compiler& c) const {
    XMMVar ret = c.newXMM();
#ifdef ASMJIT_X64
// stupid x86-64 quirk: no absolute 64-bit addressing for most instructions
    if ((sysuint_t) &this->c > (sysuint_t) 0xFFFFFFFF) {
        GPVar temp = c.newGP();
        c.mov(temp, imm((sysint_t) (void*) &this->c));
        c.movss(ret, dword_ptr(temp));
        c.unuse(temp);
    } else
#endif
        c.movss(ret, dword_ptr_abs((void*) &this->c));
    c.cvtss2sd(ret, ret); // convert float to double
    return ret;
}

XMMVar Neg::evaluate(Compiler& c) const {
    static const __v4si neg = {0, 1 << 31, 0, 1 << 31};
    XMMVar ret = a->evaluate(c);
    c.xorps(ret, dqword_ptr_abs((void*) &neg));
    return ret;
}

typedef double (*UnaryDoubleFunc)(double);
typedef double (*BinaryDoubleFunc)(double, double);
typedef double (*IntDoubleFunc)(int, double);

#define CALL_UNARY(func, arg, ret) \
    GPVar addr = c.newGP(); \
    c.mov(addr, imm((sysint_t)UnaryDoubleFunc(&func))); \
    ECall* ctx = c.call(addr); \
    c.unuse(addr); \
    ctx->setPrototype(CALL_CONV_DEFAULT, FunctionBuilder1<double, double>()); \
    XMMVar ret = c.newXMM(); \
    ctx->setArgument(0, arg); \
    ctx->setReturn(ret); \
    c.unuse(arg);

XMMVar Exp::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(std::exp, _a, ret);
    return ret;
}

XMMVar Ln::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(std::log, _a, ret);
    return ret;
}

XMMVar Sqrt::evaluate(Compiler& c) const {
    XMMVar ret = a->evaluate(c);
    c.sqrtsd(ret, ret);
    return ret;
}

XMMVar Sin::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(std::sin, _a, ret);
    return ret;
}

XMMVar Cos::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(std::cos, _a, ret);
    return ret;
}

XMMVar Tan::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(std::tan, _a, ret);
    return ret;
}

XMMVar Asin::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(std::asin, _a, ret);
    return ret;
}

XMMVar Acos::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(std::acos, _a, ret);
    return ret;
}

XMMVar Atan::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(std::atan, _a, ret);
    return ret;
}

XMMVar Gamma::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c);
    CALL_UNARY(gsl_sf_gamma, _a, ret);
    return ret;
}

XMMVar Add::evaluate(Compiler& c) const {
    XMMVar ret = a->evaluate(c),
           xb = b->evaluate(c);
    c.addsd(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Sub::evaluate(Compiler& c) const {
    XMMVar ret = a->evaluate(c),
           xb = b->evaluate(c);
    c.subsd(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Mul::evaluate(Compiler& c) const {
    XMMVar ret = a->evaluate(c),
           xb = b->evaluate(c);
    c.mulsd(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Div::evaluate(Compiler& c) const {
    XMMVar ret = a->evaluate(c),
           xb = b->evaluate(c);
    c.divsd(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar PowInt::evaluate(Compiler& c) const {
    static const double one = 1.;
    if (b == 0) {
        XMMVar ret = c.newXMM();
        c.movsd(ret, qword_ptr_abs((void*) &one));
        return ret;
    }
    XMMVar _a = a->evaluate(c),
           ret;
    int _b = b;
    bool neg = _b < 0;
    if (neg) _b = -_b;
    bool started = false;
    while (_b > 0) {
        if (_b & 1) {
            if (started) {
                c.mulsd(ret, _a);
            } else {
                if (_b == 1) {
                    ret = _a;
                } else {
                    ret = c.newXMM();
                    c.movsd(ret, _a);
                    started = true;
                }
            }
        }
        _b >>= 1;
        if (_b > 0) c.mulsd(_a, _a);
    }
    if (started) c.unuse(_a);
    if (neg) {
        XMMVar temp = c.newXMM();
        c.movsd(temp, qword_ptr_abs((void*) &one));
        c.divsd(temp, ret);
        c.unuse(ret);
        ret = temp;
    }
    return ret;
}

XMMVar Polynomial::evaluate(Compiler& c) const {
    if (!left.get()) return right->evaluate(c);
    XMMVar _l = left->evaluate(c);
    XMMVar _c = c.newXMM();
    switch (var.id->type) {
        case Variable::Id::Constant:
        case Variable::Id::Vector:
            c.movss(_c, dword_ptr_abs((void*) var.id->p));
            c.cvtss2sd(_c, _c);
            break;
        default: throw var;
    }
    c.mulsd(_l, _c);
    c.unuse(_c);
    XMMVar _r = right->evaluate(c);
    c.addsd(_l, _r);
    c.unuse(_r);
    return _l;
}

XMMVar Pow::evaluate(Compiler& c) const {
    XMMVar _a = a->evaluate(c),
           _b = b->evaluate(c);
    GPVar addr = c.newGP();
    c.mov(addr, imm((sysint_t)BinaryDoubleFunc(&std::pow)));
    ECall* ctx = c.call(addr);
    c.unuse(addr);
    ctx->setPrototype(CALL_CONV_DEFAULT, FunctionBuilder2<double, double, double>());
    XMMVar ret = c.newXMM();
    ctx->setArgument(0, _a);
    ctx->setArgument(1, _b);
    ctx->setReturn(ret);
    c.unuse(_a);
    c.unuse(_b);
    return _b;
}

XMMVar PolyGamma::evaluate(Compiler& c) const {
    XMMVar ret = a->evaluate(c);
    GPVar addr = c.newGP();
    c.mov(addr, imm((sysint_t)IntDoubleFunc(&gsl_sf_psi_n)));
    ECall* ctx = c.call(addr);
    c.unuse(addr);
    ctx->setPrototype(CALL_CONV_DEFAULT, FunctionBuilder2<double, int, double>());
    ctx->setArgument(0, imm(b));
    ctx->setArgument(1, ret);
    ctx->setReturn(ret);
    return ret;
}

#if 0
VectorFunc Expression::evaluatorVector() const {
    Compiler c;
    c.setLogger(new FileLogger(stdout));
    c.newFunction(CALL_CONV_DEFAULT, FunctionBuilder1<float, int>()); // HACK: returning a "float" but really a v4sf
    c.ret(evaluateVector(c, c.argGP(0)));
    c.endFunction();
    return (VectorFunc)c.make();
}

XMMVar Variable::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = c.newXMM();
    switch (id->type) {
        case Id::Constant:
#ifdef ASMJIT_X64
// stupid x86-64 quirk: no absolute 64-bit addressing for most instructions
            if ((sysuint_t) id->p > (sysuint_t) 0xFFFFFFFF) {
                GPVar temp = c.newGP();
                c.mov(temp, imm((sysint_t) (void*) id->p));
                c.movss(ret, dword_ptr(temp));
                c.unuse(temp);
            } else
#endif
                c.movss(ret, dword_ptr_abs((void*) id->p));
            c.shufps(ret, ret, imm(0x00)); // replicate least significant data element
            return ret;
        case Id::Vector:
#ifdef ASMJIT_X64
// stupid x86-64 quirk: no absolute 64-bit addressing for most instructions
            if ((sysuint_t) id->p > (sysuint_t) 0xFFFFFFFF) {
                GPVar temp = c.newGP();
                c.mov(temp, imm((sysint_t) (void*) id->p));
                c.movss(ret, dword_ptr(temp, i, 2));
                c.unuse(temp);
            } else
#endif
                c.movaps(ret, dqword_ptr_abs((void*) id->p, i, 2));
            return ret;
        default:
            throw this;
    }
}

XMMVar Constant::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = c.newXMM();
#ifdef ASMJIT_X64
// stupid x86-64 quirk: no absolute 64-bit addressing for most instructions
    if ((sysuint_t) &this->c > (sysuint_t) 0xFFFFFFFF) {
        GPVar temp = c.newGP();
        c.mov(temp, imm((sysint_t) (void*) &this->c));
        c.movss(ret, dword_ptr(temp));
        c.unuse(temp);
    } else
#endif
        c.movss(ret, dword_ptr_abs((void*) &this->c));
    c.shufps(ret, ret, imm(0x00)); // replicate least significant data element
    return ret;
}

XMMVar Neg::evaluateVector(Compiler& c, const GPVar& i) const {
    static const __v4si neg = {1 << 31, 1 << 31, 1 << 31, 1 << 31};
    XMMVar ret = a->evaluateVector(c, i);
    c.xorps(ret, dqword_ptr_abs((void*) &neg));
    return ret;
}

#define CALL_SSE(func, arg, ret) \
    GPVar addr = c.newGP(); \
    c.mov(addr, imm((sysint_t)(&func))); \
    ECall* ctx = c.call(addr); \
    c.unuse(addr); \
    ctx->setPrototype(CALL_CONV_DEFAULT, FunctionBuilder1<float, float>()); /* HACK */ \
    XMMVar ret = c.newXMM(); \
    ctx->setArgument(0, arg); \
    ctx->setReturn(ret); \
    c.unuse(arg);

XMMVar Exp::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar _a = a->evaluateVector(c, i);
    CALL_SSE(exp_ps, _a, ret);
    return ret;
}

XMMVar Ln::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar _a = a->evaluateVector(c, i);
    CALL_SSE(log_ps, _a, ret);
    return ret;
}

XMMVar Sqrt::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i);
    c.sqrtps(ret, ret);
    return ret;
}

XMMVar Sin::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar _a = a->evaluateVector(c, i);
    CALL_SSE(sin_ps, _a, ret);
    return ret;
}

XMMVar Cos::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar _a = a->evaluateVector(c, i);
    CALL_SSE(cos_ps, _a, ret);
    return ret;
}

XMMVar Tan::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar _a = a->evaluateVector(c, i);
    CALL_SSE(tan_ps, _a, ret);
    return ret;
}

#define CALL_UNARY_VECTOR(func, arg, ret)

XMMVar Asin::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Acos::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Atan::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Gamma::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar Add::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i),
           xb = b->evaluateVector(c, i);
    c.addps(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Sub::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i),
           xb = b->evaluateVector(c, i);
    c.subps(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Mul::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i),
           xb = b->evaluateVector(c, i);
    c.mulps(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar Div::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i),
           xb = b->evaluateVector(c, i);
    c.divps(ret, xb);
    c.unuse(xb);
    return ret;
}

XMMVar PowInt::evaluateVector(Compiler& c, const GPVar& i) const {
    if (b == 0) {
        XMMVar ret = c.newXMM();;
        static const __v4sf ones = {1.f, 1.f, 1.f, 1.f};
        c.movaps(ret, dqword_ptr_abs((void*) &ones));
        return ret;
    }
    XMMVar _a = a->evaluateVector(c, i),
           ret;
    int _b = b;
    bool neg = _b < 0;
    if (neg) _b = -_b;
    bool started = false;
    while (_b > 0) {
        if (_b & 1) {
            if (started) {
                c.mulps(ret, _a);
            } else {
                if (_b == 1) {
                    ret = _a;
                } else {
                    ret = c.newXMM();
                    c.movaps(ret, _a);
                    started = true;
                }
            }
        }
        _b >>= 1;
        if (_b > 0) c.mulps(_a, _a);
    }
    if (started) c.unuse(_a);
    if (neg) c.rcpps(ret, ret);
    return ret;
}

XMMVar Polynomial::evaluateVector(Compiler& c, const GPVar& i) const {
    if (!left.get()) return right->evaluateVector(c, i);
    XMMVar _l = left->evaluateVector(c, i);
    XMMVar _c = c.newXMM();
    switch (var.id->type) {
        case Variable::Id::Constant: {
            c.movss(_c, dword_ptr_abs((void*) var.id->p));
            c.shufps(_c, _c, imm(0x00));
            break;
        }
        case Variable::Id::Vector: {
            c.movaps(_c, dqword_ptr_abs((void*) var.id->p, i, 2));
            break;
        }
        default: throw var;
    }
    c.mulps(_l, _c);
    c.unuse(_c);
    XMMVar _r = right->evaluateVector(c, i);
    c.addps(_l, _r);
    c.unuse(_r);
    return _l;
}

XMMVar Pow::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i);
    throw "Not implemented";
    return ret;
}

XMMVar PolyGamma::evaluateVector(Compiler& c, const GPVar& i) const {
    XMMVar ret = a->evaluateVector(c, i);
    throw "Not implemented";
    return ret;
}
#endif
