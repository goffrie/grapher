#ifndef _EXPRESSION_H_
#define _EXPRESSION_H_

#include <boost/lexical_cast.hpp>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <memory>
#include <unordered_map>
#include <QMetaType>

typedef double Number;
typedef Number* Vector;
typedef Number* __restrict__ VectorR;

struct Variable;
namespace std {
    template<> class hash<Variable> : public std::unary_function<Variable, size_t> {
    public:
        std::size_t operator()(const Variable& a) const;
    };
}
namespace Precedence {
    enum {
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
    EPtr ecopy() const { return EPtr(copy()); }
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
    std::string toString(int prec = -1) const { return boost::lexical_cast<std::string>(c); }
};
struct External : public Expression {
    Vector c;
    External(Vector _c) : c(_c) { }
    static std::unique_ptr<External> create(Vector _c) { return std::unique_ptr<External>(new External(_c)); }
    Number evaluate() const { return *c; }
    Vector evaluateVector(std::size_t size) const;
    EPtr substitute(const Subst&) const { return EPtr(copy()); }
    EPtr derivative(const Variable&) const { return EPtr(new Constant(0)); }
	void variables(std::set<Variable>&) const { }
    External* copy() const { return new External(c); }
    std::string toString(int prec = -1) const { return std::string("External(") + boost::lexical_cast<std::string>(c) + ")"; }
};
struct Variable : public Expression {
	typedef std::string Id;
	std::shared_ptr<Id> id;
    Variable() : id(new Id) { }
    Variable(Id _id) : id(new Id(_id)) { }
	Variable(const Variable& b) : id(b.id) { }
	Variable(const std::shared_ptr<Id>& b) : id(b) { }
    Number evaluate() const { throw this; }
    Vector evaluateVector(std::size_t size) const { throw this; }
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
	bool operator==(const Variable& b) const { return id == b.id; }
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
    EPtr simplify() const;
	void variables(std::set<Variable>& vars) const { a->variables(vars); }
	virtual UnaryOp* construct(EPtr _a) const = 0;
    EPtr substitute(const Subst& s) const { return EPtr(construct(a->substitute(s))); }
    UnaryOp* copy() const { return construct(EPtr(a->copy())); }
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

UNARY_FUNCTION(Neg, -, -, EPtr simplify() const; \
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
};

#define BINARY_OP(Name, op, ladj, radj, func) \
struct Name : public BinaryOp { \
    Name(EPtr _a, EPtr _b) : BinaryOp(std::move(_a), std::move(_b)) { } \
    Number evaluate() const { return func(a->evaluate(), b->evaluate()); } \
    Number evaluate(Number _a, Number _b) const { return func(_a, _b); } \
    Vector evaluateVector(std::size_t size) const; \
    EPtr derivative(const Variable& var) const; \
    EPtr simplify() const; \
    Name* construct(EPtr _a, EPtr _b) const { return new Name(std::move(_a), std::move(_b)); } \
    static std::unique_ptr<Name> create(EPtr _a, EPtr _b) { return std::unique_ptr<Name>(new Name(std::move(_a), std::move(_b))); } \
    std::string toString(int prec = -1) const { return wrap(prec, Precedence::Name, a->toString(Precedence::Name+ladj) + " " #op " " + b->toString(Precedence::Name+radj)); } \
};

#define SIMPLE_OP(Name, op, ladj, radj) BINARY_OP(Name, op, ladj, radj, [](Number a, Number b) -> Number { return a op b; })

SIMPLE_OP(Add, +, 0, 0)
SIMPLE_OP(Sub, -, 0, 1)
SIMPLE_OP(Mul, *, 0, 0)
SIMPLE_OP(Div, /, 0, 1)
BINARY_OP(Pow, ^, 1, 0, std::pow)
/*struct Add : public BinaryOp {
    Add(EPtr _a, EPtr _b) : BinaryOp(std::move(_a), std::move(_b)) { }
    Number evaluate() const { return a->evaluate() + b->evaluate(); }
    Number evaluate(Number _a, Number _b) const { return _a + _b; }
    Vector evaluateVector(std::size_t size) const;
    EPtr derivative(const Variable& var) const;
    EPtr simplify() const;
    Add* construct(EPtr _a, EPtr _b) const { return new Add(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") + (" + b->toString() + ")"; }
};
struct Sub : public BinaryOp {
    Sub(EPtr _a, EPtr _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return a->evaluate() - b->evaluate(); }
    Number evaluate(Number _a, Number _b) const { return _a - _b; }
    Vector evaluateVector(std::size_t size) const;
    EPtr derivative(const Variable& var) const;
    EPtr simplify() const;
    Sub* construct(EPtr _a, EPtr _b) const { return new Sub(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") - (" + b->toString() + ")"; }
};
struct Mul : public BinaryOp {
    Mul(EPtr _a, EPtr _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return a->evaluate() * b->evaluate(); }
    Number evaluate(Number _a, Number _b) const { return _a * _b; }
    Vector evaluateVector(std::size_t size) const;
    EPtr derivative(const Variable& var) const;
    EPtr simplify() const;
    Mul* construct(EPtr _a, EPtr _b) const { return new Mul(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") * (" + b->toString() + ")"; }
};
struct Div : public BinaryOp {
    Div(EPtr _a, EPtr _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return a->evaluate() / b->evaluate(); }
    Number evaluate(Number _a, Number _b) const { return _a / _b; }
    Vector evaluateVector(std::size_t size) const;
    EPtr derivative(const Variable& var) const;
    EPtr simplify() const;
    Div* construct(EPtr _a, EPtr _b) const { return new Div(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") / (" + b->toString() + ")"; }
};
struct Pow : public BinaryOp {
    Pow(EPtr _a, EPtr _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return std::pow(a->evaluate(), b->evaluate()); }
    Number evaluate(Number _a, Number _b) const { return std::pow(_a, _b); }
    Vector evaluateVector(std::size_t size) const;
    EPtr derivative(const Variable& var) const;
    EPtr simplify() const;
    Pow* construct(EPtr _a, EPtr _b) const { return new Pow(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") ^ (" + b->toString() + ")"; }
};*/
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
    std::string toString(int prec = -1) const { return wrap(prec, Precedence::Pow, a->toString(Precedence::Pow+1) + " ^ !" + boost::lexical_cast<std::string>(b)); }
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
    std::string toString(int prec = -1) const { return wrap(prec, Precedence::Func, std::string("psi_") + boost::lexical_cast<std::string>(b) + "(" + a->toString(-1) + ")"); }
};
struct Equation : public Thing {
    EPtr a, b;
    Equation(EPtr _a, EPtr _b) : a(std::move(_a)), b(std::move(_b)) { }
    static std::unique_ptr<Equation> create(EPtr _a, EPtr _b) { return std::unique_ptr<Equation>(new Equation(std::move(_a), std::move(_b))); }
    std::string toString(int prec = -1) const { return wrap(prec, Precedence::Eq, a->toString(Precedence::Eq) + " = " + b->toString(Precedence::Eq)); }
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

#endif
