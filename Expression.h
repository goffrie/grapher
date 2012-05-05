#ifndef _EXPRESSION_H_
#define _EXPRESSION_H_

#include <exception>

#include <memory>
#include <unordered_map>
#include <vector>
#include <set>
#include <string>
#include <cstdint>

#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <boost/config/suffix.hpp>

#include <QMetaType>

#include <malloc.h>

inline void* aligned_malloc(size_t size) {
    void* pa = malloc(size + 15 + sizeof(void *));
    if (!pa) return NULL;
    void* ptr = (void*) ( ((uintptr_t)pa + sizeof(void *) + 15) & (~15) );
    *((void **)ptr-1) = pa;
    return ptr;
}

inline void aligned_free(void* ptr) {
    if (ptr) free(*((void **)ptr-1));
}

typedef float Number;

BOOST_CONSTEXPR_OR_CONST int SSE_VECTOR_SIZE = 4;

#define VECTOR_ALLOC(num) ((Vector) aligned_malloc((((num) + SSE_VECTOR_SIZE - 1) & (~(SSE_VECTOR_SIZE - 1)))*sizeof(Number)))
#define VECTOR_FREE(ptr) aligned_free(ptr)

typedef Number __attribute__((aligned(16))) * Vector;
typedef Number __attribute__((aligned(16))) * __restrict VectorR;

struct VectorDeleter {
    void operator()(Vector ptr) const throw() {
        VECTOR_FREE(ptr);
    }
};

typedef std::unique_ptr<Number[], VectorDeleter> UVector;

static_assert(sizeof(UVector) == sizeof(Vector), "UVector is not the right size!");

struct Variable;
namespace std {
    template<> class hash<Variable> : public std::unary_function<Variable, size_t> {
    public:
        std::size_t operator()(const Variable& a) const;
    };
}
namespace Precedence {
    enum {
        Const = 1000,
        Func = 100,
        Pow = 10,
        Div = 7,
        Mul = 7,
        Neg = 5,
        Sub = 3,
        Add = 3,
        Eq = 0
    };
}

struct Thing {
    virtual ~Thing() { };
    virtual std::string toString(int precedence = -1) const = 0;
};
struct Expression;
struct Polynomial;
typedef std::unique_ptr<Expression> EPtr;
struct Expression : public Thing {
    typedef std::unordered_map<Variable, Expression*> Subst;
    virtual Number evaluate() const = 0;
    virtual Vector evaluateVector(std::size_t size) const = 0;
    virtual EPtr substitute(const Subst& s) const = 0;
    virtual EPtr derivative(const Variable& var) const = 0;
    virtual EPtr simplify() const { return EPtr(copy()); }
    virtual void variables(std::set<Variable>& out) const = 0;
    virtual Expression* copy() const = 0;
    virtual bool polynomial(const Variable&) const { return false; }
    EPtr ecopy() const { return EPtr(copy()); }
    virtual EPtr expand() const { return ecopy(); }
    virtual std::unique_ptr<Polynomial> facsum(const Variable& var) const;
};
struct Constant : public Expression {
    Number c;
    Constant(Number _c) : c(_c) { }
    static std::unique_ptr<Constant> create(Number _c) { return std::unique_ptr<Constant>(new Constant(_c)); }
    Number evaluate() const { return c; }
    Vector evaluateVector(std::size_t size) const;
    EPtr substitute(const Subst&) const { return EPtr(copy()); }
    EPtr derivative(const Variable&) const { return EPtr(new Constant(0)); }
    void variables(std::set<Variable>&) const { }
    Constant* copy() const { return new Constant(c); }
    bool polynomial(const Variable&) const { return true; }
    std::string toString(int prec = -1) const;
};
struct Variable : public Expression {
    struct Id {
        std::string s;
        enum Type { None, Constant, Vector } type;
        ::Vector p;
        Id() : type(None) { }
        Id(const char* _s) : s(_s), type(None) { }
        Id(std::string _s) : s(_s), type(None) { }
        Id(std::string _s, Type t, ::Vector _p = NULL) : s(_s), type(t), p(_p) { }
        operator std::string() const { return s; }
    };
    std::shared_ptr<Id> id;
    Variable() : id(new Id) { }
    Variable(Id _id) : id(new Id(_id)) { }
    Variable(const Variable& b) : id(b.id) { }
    Variable(const std::shared_ptr<Id>& b) : id(b) { }
    Variable& operator=(const Variable& b) = default;
    static std::unique_ptr<Variable> create() { return std::unique_ptr<Variable>(new Variable); }
    static std::unique_ptr<Variable> create(const Variable& b) { return std::unique_ptr<Variable>(new Variable(b)); }
    Number evaluate() const;
    Vector evaluateVector(std::size_t size) const;
    EPtr substitute(const Subst& s) const {
        Subst::const_iterator p = s.find(*this);
        if (p == s.end()) return EPtr(copy());
        return EPtr(p->second->copy());
    }
    EPtr derivative(const Variable& var) const {
        return EPtr(new Constant((var == *this) ? 1.0 : 0.0));
    }
    void variables(std::set<Variable>& vars) const { vars.insert(*this); }
    Variable* copy() const { return new Variable(id); }
    bool polynomial(const Variable&) const { return true; }
    std::unique_ptr<Polynomial> facsum(const Variable& var) const;
    bool operator==(const Variable& b) const { return id == b.id; }
    bool operator!=(const Variable& b) const { return id != b.id; }
    friend bool operator<(const Variable& a, const Variable& b) { return a.id < b.id; }
    std::string toString(int prec = -1) const { return *id; }
};
Q_DECLARE_METATYPE(Variable);
namespace std {
    inline std::size_t hash<Variable>::operator()(const Variable& a) const {
        return std::hash<Variable::Id*>()(a.id.get());
    }
}
struct UnaryOp : public Expression {
    EPtr a;
    UnaryOp(EPtr _a) : a(std::move(_a)) { }
    virtual Number evaluate(Number _a) const = 0;
    virtual Number evaluate() const override = 0;
    EPtr simplify() const;
    void variables(std::set<Variable>& vars) const { a->variables(vars); }
    virtual UnaryOp* construct(EPtr _a) const = 0;
    EPtr substitute(const Subst& s) const { return EPtr(construct(a->substitute(s))); }
    UnaryOp* copy() const { return construct(EPtr(a->copy())); }
    EPtr expand() const { return EPtr(construct(a->expand())); }
};

#define wrap(P, Q, S) (((P) > (Q)) ? std::string("(") : std::string()) + (S) + (((P) > (Q)) ? std::string(")") : std::string())

#define UNARY_FUNCTION(Name, func, sfunc, extra) \
struct Name : public UnaryOp { \
    Name(EPtr _a) : UnaryOp(std::move(_a)) { } \
    Number evaluate(Number _a) const; \
    Number evaluate() const; \
    Vector evaluateVector(std::size_t size) const; \
    EPtr derivative(const Variable& var) const; \
    Name* construct(EPtr _a) const { return new Name(std::move(_a)); } \
    static std::unique_ptr<Name> create(EPtr _a) { return std::unique_ptr<Name>(new Name(std::move(_a))); } \
    extra \
};
#define FUNCTION_PRINTER(func) std::string toString(int prec = -1) const { return wrap(prec, Precedence::Func, std::string(#func "(") + a->toString(-1) + ")"); }
#define SIMPLE_UNARY_FUNCTION(Name, func) UNARY_FUNCTION(Name, func, std::func, FUNCTION_PRINTER(func) )

UNARY_FUNCTION(Neg, -, -, bool polynomial(const Variable& var) const { return a->polynomial(var); } \
EPtr simplify() const; \
EPtr expand() const; \
std::unique_ptr<Polynomial> facsum(const Variable&) const; \
std::string toString(int prec = -1) const { return std::string("-") + a->toString(Precedence::Neg); })
SIMPLE_UNARY_FUNCTION(Exp, exp)
SIMPLE_UNARY_FUNCTION(Ln, log)
SIMPLE_UNARY_FUNCTION(Sqrt, sqrt)
SIMPLE_UNARY_FUNCTION(Sin, sin)
SIMPLE_UNARY_FUNCTION(Cos, cos)
SIMPLE_UNARY_FUNCTION(Tan, tan)
SIMPLE_UNARY_FUNCTION(Asin, asin)
SIMPLE_UNARY_FUNCTION(Acos, acos)
SIMPLE_UNARY_FUNCTION(Atan, atan)
UNARY_FUNCTION(Gamma, gamma, gsl_sf_gamma, FUNCTION_PRINTER(gamma) )

#undef UNARY_FUNCTION

struct BinaryOp : public Expression {
    EPtr a, b;
    BinaryOp(EPtr _a, EPtr _b) : a(std::move(_a)), b(std::move(_b)) { }
    void variables(std::set<Variable>& vars) const { a->variables(vars); b->variables(vars); }
    virtual BinaryOp* construct(EPtr _a, EPtr _b) const = 0;
    EPtr substitute(const Subst& s) const { return EPtr(construct(EPtr(a->substitute(s)), EPtr(b->substitute(s)))); }
    BinaryOp* copy() const { return construct(EPtr(a->copy()), EPtr(b->copy())); }
    EPtr expand() const { return EPtr(construct(a->expand(), b->expand())); }
};

#define BINARY_OP(Name, op, ladj, radj, expr, extra) \
struct Name : public BinaryOp { \
    Name(EPtr _a, EPtr _b) : BinaryOp(std::move(_a), std::move(_b)) { } \
    Number evaluate() const { return evaluate(a->evaluate(), b->evaluate()); } \
    Number evaluate(Number _1, Number _2) const { return expr; } \
    Vector evaluateVector(std::size_t size) const; \
    EPtr derivative(const Variable& var) const; \
    EPtr simplify() const; \
    Name* construct(EPtr _a, EPtr _b) const { return new Name(std::move(_a), std::move(_b)); } \
    static std::unique_ptr<Name> create(EPtr _a, EPtr _b) { return std::unique_ptr<Name>(new Name(std::move(_a), std::move(_b))); } \
    std::string toString(int prec = -1) const { return wrap(prec, Precedence::Name, a->toString(Precedence::Name+ladj) + " " #op " " + b->toString(Precedence::Name+radj)); } \
    extra \
};

#define SIMPLE_OP(Name, op, ladj, radj, extra) \
    BINARY_OP(Name, op, ladj, radj, \
        _1 op _2, \
        bool polynomial(const Variable& var) const { return a->polynomial(var) && b->polynomial(var); } \
        EPtr expand() const; \
        extra)

SIMPLE_OP(Add, +, 0, 0, std::unique_ptr<Polynomial> facsum(const Variable&) const;)
SIMPLE_OP(Sub, -, 0, 1, )
SIMPLE_OP(Mul, *, 0, 0, std::unique_ptr<Polynomial> facsum(const Variable&) const;)
BINARY_OP(Div, /, 0, 1, _1 / _2, \
        bool polynomial(const Variable& var) const; \
        EPtr expand() const; \
        std::unique_ptr<Polynomial> facsum(const Variable&) const;)
BINARY_OP(Pow, ^, 1, 0, std::pow(_1, _2), bool polynomial(const Variable& var) const;)

struct PowInt : public UnaryOp {
    int b;
    PowInt(EPtr _a, int _b) : UnaryOp(std::move(_a)), b(_b) { }
    static std::unique_ptr<PowInt> create(EPtr a, int b) { return std::unique_ptr<PowInt>(new PowInt(std::move(a), b)); }
    Number evaluate() const;
    Number evaluate(Number _a) const;
    Vector evaluateVector(std::size_t size) const;
    EPtr derivative(const Variable& var) const;
    EPtr simplify() const;
    PowInt* construct(EPtr _a) const { return new PowInt(std::move(_a), b); }
    bool polynomial(const Variable& var) const { return a->polynomial(var) && b >= 0; }
    EPtr expand() const;
    std::unique_ptr<Polynomial> facsum(const Variable& var) const;
    std::string toString(int prec = -1) const;
};

// left * var + right
struct Polynomial : public Expression {
    std::unique_ptr<Polynomial> left;
    Variable var;
    EPtr right;
    Polynomial(std::unique_ptr<Polynomial> l, const Variable& v, EPtr r) : left(std::move(l)), var(v), right(std::move(r)) {
        Q_ASSERT(left.get() == NULL || left->var == var);
    }
    Polynomial(const Variable& v, EPtr r) : var(v), right(std::move(r)) { }
    Number evaluate() const;
    Vector evaluateVector(std::size_t size) const;
    EPtr substitute(const Subst& s) const;
    EPtr derivative(const Variable& var) const;
    EPtr simplify() const;
    void variables(std::set<Variable>& out) const;
    Polynomial* copy() const;
    bool polynomial(const Variable&) const;
    std::unique_ptr<Polynomial> facsum(const Variable& var) const;
    int degree() const;
    static EPtr create(EPtr l, const Variable& v, EPtr r);
    static std::unique_ptr<Polynomial> create(std::unique_ptr<Polynomial> l, const Variable& v, EPtr r);
    static std::unique_ptr<Polynomial> create(const Variable& v, EPtr r);
    std::string toString(int prec = -1) const;
};

struct PolyGamma : public UnaryOp {
    int b;
    PolyGamma(EPtr _a, int _b) : UnaryOp(std::move(_a)), b(_b) { }
    static std::unique_ptr<PolyGamma> create(EPtr a, int b) { return std::unique_ptr<PolyGamma>(new PolyGamma(std::move(a), b)); }
    Number evaluate() const;
    Number evaluate(Number _a) const;
    Vector evaluateVector(std::size_t size) const;
    EPtr derivative(const Variable& var) const;
    PolyGamma* construct(EPtr _a) const { return new PolyGamma(std::move(_a), b); }
    std::string toString(int prec = -1) const;
};
struct Equation : public Thing {
    EPtr a, b;
    Equation(EPtr _a, EPtr _b) : a(std::move(_a)), b(std::move(_b)) { }
    static std::unique_ptr<Equation> create(EPtr _a, EPtr _b) { return std::unique_ptr<Equation>(new Equation(std::move(_a), std::move(_b))); }
    std::string toString(int prec = -1) const { return wrap(prec, Precedence::Eq, a->toString(Precedence::Eq) + " = " + b->toString(Precedence::Eq)); }
};
struct Inequality : public Thing {
    EPtr a, b;
    enum Type {
        LT,
        GT,
        LTE,
        GTE
    };
    Type type;
    struct UnknownSignException : public std::exception { const char* what() const throw() { return "Unknown inequality sign"; } };
    static std::string sign(Type t) {
        switch (t) {
            case LT: return "<";
            case GT: return ">";
            case LTE: return "<=";
            case GTE: return ">=";
        }
        throw UnknownSignException();
    }
    static Type fromSign(std::string t) {
        if (t == "<") return LT;
        if (t == ">") return GT;
        if (t == "<=") return LTE;
        if (t == ">=") return GTE;
        throw UnknownSignException();
    }
    Inequality(EPtr _a, EPtr _b, Type _t) : a(std::move(_a)), b(std::move(_b)), type(_t) { }
    static std::unique_ptr<Inequality> create(EPtr _a, EPtr _b, Type _t) { return std::unique_ptr<Inequality>(new Inequality(std::move(_a), std::move(_b), _t)); }
    std::unique_ptr<Inequality> substitute(const Expression::Subst& s) const { return create(a->substitute(s), b->substitute(s), type); }
    std::unique_ptr<Inequality> simplify() const;
    std::string toString(int prec = -1) const { return wrap(prec, Precedence::Eq, a->toString(Precedence::Eq) + " " + sign(type) + " " + b->toString(Precedence::Eq)); }
    Vector evaluateVector(std::size_t size) const;
};
/*
struct Function : public Thing {
    struct ArgumentMismatchException { };
    struct UnboundVariableException { };
    std::vector<Variable> bound_args;
    EPtr body;
    Function() : body(NULL) { }
    Function(std::vector<Variable> _args, EPtr _body) : bound_args(_args), body(_body) { check(); }
    ~Function() { delete body; }
    void check() {
        std::set<Variable> vars;
        body->variables(vars);
        if (vars.size() != bound_args.size()) throw UnboundVariableException();
        for (std::size_t i = 0; i < bound_args.size(); ++i) {
            if (vars.find(bound_args[i]) == vars.end()) throw UnboundVariableException();
        }
    }
    EPtr call(const std::vector<EPtr>& args) {
        if (args.size() != bound_args.size()) throw ArgumentMismatchException();
        Expression::Subst s;
        for (std::size_t i = 0; i < bound_args.size(); ++i) {
            s.insert(std::make_pair(bound_args[i], args[i]));
        }
        return body->substitute(s);
    }
};
*/

inline EPtr operator-(EPtr a) {
    return Neg::create(std::move(a));
}
inline EPtr operator+(EPtr a, EPtr b) {
    return Add::create(std::move(a), std::move(b));
}
inline EPtr operator-(EPtr a, EPtr b) {
    return Sub::create(std::move(a), std::move(b));
}
inline EPtr operator*(EPtr a, EPtr b) {
    return Mul::create(std::move(a), std::move(b));
}
inline EPtr operator/(EPtr a, EPtr b) {
    return Div::create(std::move(a), std::move(b));
}

#endif
