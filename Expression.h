#ifndef _EXPRESSION_H_
#define _EXPRESSION_H_

#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <set>
#include <string>
#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

// #define USE_LLVM

#ifdef USE_LLVM
#include <llvm/DerivedTypes.h>
#include <llvm/LLVMContext.h>
#include <llvm/Module.h>
#include <llvm/Analysis/Verifier.h>
#include <llvm/Support/IRBuilder.h>
#endif

typedef float Number;
// typedef void (*Evaluator)(Number*);
typedef Number* Vector;
typedef Number* __restrict__ VectorR;

struct Variable;
struct Thing {
    virtual ~Thing() { };
    virtual std::string toString() const = 0;
};
struct Expression : public Thing {
    typedef boost::unordered_map<Variable, Expression*> Subst;
    virtual Number evaluate() const = 0;
    virtual Vector evaluateVector(std::size_t size) const = 0;
#ifdef USE_LLVM
    // the resulting function is bound to the lifetime of this object
    Evaluator evaluator() const {
        AsmJit::Compiler c;
        c.newFunction(AsmJit::CALL_CONV_DEFAULT, AsmJit::FunctionBuilder1<AsmJit::Void, Number*>());
        c.getFunction()->setHint(AsmJit::FUNCTION_HINT_NAKED, true);
        AsmJit::XMMVar retvar(evaluate(c));
        c.movss(AsmJit::Mem(c.argGP(0), 0), retvar);
        c.endFunction();
        return AsmJit::function_cast<Evaluator>(c.make());
    }
    virtual llvm::Value* evaluate(llvm::IRBuilder<>& c) const = 0;
#endif
    virtual Expression* substitute(const Subst& s) const = 0;
    virtual Expression* derivative(const Variable& var) const = 0;
    virtual Expression* simplify() const { return copy(); }
	virtual void variables(std::set<Variable>& out) const = 0;
    virtual Expression* copy() const = 0;
};
struct Constant : public Expression {
    Number c;
    Constant(Number _c) : c(_c) { }
    Number evaluate() const { return c; }
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& cmp) const {
        return llvm::ConstantFP::get(llvm::getGlobalContext(), c);
    }
#endif
    Expression* substitute(const Subst&) const { return copy(); }
    Constant* derivative(const Variable& /*var*/) const { return new Constant(0); }
	void variables(std::set<Variable>&) const { }
    Constant* copy() const { return new Constant(c); }
    std::string toString() const { return boost::lexical_cast<std::string>(c); }
};
struct External : public Expression {
    Number* __restrict__ c;
    External(Number* _c) : c(_c) { }
    Number evaluate() const { return *c; }
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& cmp) const {
        return cmp.CreateLoad(llvm::GlobalVariable());
    }
#endif
    Expression* substitute(const Subst&) const { return copy(); }
    Constant* derivative(const Variable& /*var*/) const { return new Constant(0); }
	void variables(std::set<Variable>&) const { }
    External* copy() const { return new External(c); }
    std::string toString() const { return std::string("External(") + boost::lexical_cast<std::string>(c) + ")"; }
};
struct Variable : public Expression {
	typedef std::string Id;
	boost::shared_ptr<Id> id;
    Variable() : id(new Id) { }
    Variable(Id _id) : id(new Id(_id)) { }
	Variable(const Variable& b) : id(b.id) { }
	Variable(const boost::shared_ptr<Id>& b) : id(b) { }
    Number evaluate() const { throw this; }
    Vector evaluateVector(std::size_t size) const { throw this; }
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>&) const { throw this; }
#endif
    Expression* substitute(const Subst& s) const {
        Subst::const_iterator p = s.find(*this);
        if (p == s.end()) return copy();
        return p->second->copy();
    }
    Constant* derivative(const Variable& var) const {
        return new Constant((var == *this) ? 1.0 : 0.0);
    }
	void variables(std::set<Variable>& vars) const { vars.insert(*this); }
    Variable* copy() const { return new Variable(id); }
	bool operator==(const Variable& b) const { return id == b.id; }
	friend std::size_t hash_value(const Variable& a) { return boost::hash<boost::shared_ptr<Id> >()(a.id); }
	friend bool operator<(const Variable& a, const Variable& b) { return a.id < b.id; }
    std::string toString() const { return *id + "#" + boost::lexical_cast<std::string>(id.get()); }
};
struct UnaryOp : public Expression {
    Expression* a;
    UnaryOp(Expression* _a) : a(_a) { }
    ~UnaryOp() { delete a; }
    virtual Number evaluate(Number _a) const = 0;
    Expression* simplify() const;
	void variables(std::set<Variable>& vars) const { a->variables(vars); }
	virtual UnaryOp* construct(Expression* _a) const = 0;
    UnaryOp* copy() const { return construct(a->copy()); }
};
struct Neg : public UnaryOp {
    Neg(Expression* _a) : UnaryOp(_a) { }
    Number evaluate(Number _a) const { return -_a; }
    Number evaluate() const { return -a->evaluate(); }
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& c) const {
        return c.CreateFNeg(a->evaluate(c));
    }
#endif
    Expression* substitute(const Subst& s) const { return new Neg(a->substitute(s)); }
    Expression* derivative(const Variable& var) const { return new Neg(a->derivative(var)); }
    Expression* simplify() const;
    Neg* copy() const { return new Neg(a->copy()); }
    Neg* construct(Expression* _a) const { return new Neg(_a); }
    std::string toString() const { return std::string("-(") + a->toString() + ")"; }
};
#define UNARY_FUNCTION(Name, func, sfunc) \
struct Name : public UnaryOp { \
    Name(Expression* _a) : UnaryOp(_a) { } \
    Number evaluate(Number _a) const { return sfunc(_a); } \
    Number evaluate() const { return sfunc(a->evaluate()); } \
    Vector evaluateVector(std::size_t size) const { \
        Vector _a = a->evaluateVector(size); \
        for (std::size_t i = 0; i < size; ++i) _a[i] = sfunc(_a[i]); \
        return _a; \
    } \
    Expression* substitute(const Subst& s) const { return new Name(a->substitute(s)); } \
    Expression* derivative(const Variable& var) const; \
    Name* construct(Expression* _a) const { return new Name(_a); } \
    std::string toString() const { return std::string(#func "(") + a->toString() + ")"; } \
};

/*
     llvm::Value* evaluate(llvm::IRBuilder<>& c) const { \
        AsmJit::XMMVar ret(c.newXMM()); \
        Number (*func)(Number) = &std::func; \
        AsmJit::ECall* ctx = c.call((void*)func); \
        ctx->setPrototype(AsmJit::CALL_CONV_DEFAULT, AsmJit::FunctionBuilder1<Number, Number>()); \
        ctx->setArgument(0, a->evaluate(c)); \
        ctx->setReturn(ret); \
        return ret; \
    } \
*/

// Number gamma(Number a);

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

struct BinaryOp : public Expression {
    Expression* a;
    Expression* b;
    BinaryOp(Expression* _a, Expression* _b) : a(_a), b(_b) { }
    ~BinaryOp() { delete a; delete b; }
	void variables(std::set<Variable>& vars) const { a->variables(vars); b->variables(vars); }
	virtual BinaryOp* construct(Expression* _a, Expression* _b) const = 0;
    BinaryOp* copy() const { return construct(a->copy(), b->copy()); }
};
struct Add : public BinaryOp {
    Add(Expression* _a, Expression* _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return a->evaluate() + b->evaluate(); }
    Number evaluate(Number _a, Number _b) const { return _a + _b; }
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& c) const { return c.CreateFAdd(a->evaluate(c), b->evaluate(c)); }
#endif
    Expression* substitute(const Subst& s) const { return new Add(a->substitute(s), b->substitute(s)); }
    Expression* derivative(const Variable& var) const;
    Expression* simplify() const;
    Add* construct(Expression* _a, Expression* _b) const { return new Add(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") + (" + b->toString() + ")"; }
};
struct Sub : public BinaryOp {
    Sub(Expression* _a, Expression* _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return a->evaluate() - b->evaluate(); }
    Number evaluate(Number _a, Number _b) const { return _a - _b; }
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& c) const { return c.CreateFSub(a->evaluate(c), b->evaluate(c)); }
#endif
    Expression* substitute(const Subst& s) const { return new Sub(a->substitute(s), b->substitute(s)); }
    Expression* derivative(const Variable& var) const;
    Expression* simplify() const;
    Sub* construct(Expression* _a, Expression* _b) const { return new Sub(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") - (" + b->toString() + ")"; }
};
struct Mul : public BinaryOp {
    Mul(Expression* _a, Expression* _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return a->evaluate() * b->evaluate(); }
    Number evaluate(Number _a, Number _b) const { return _a * _b; }
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& c) const { return c.CreateFMul(a->evaluate(c), b->evaluate(c)); }
#endif
    Expression* substitute(const Subst& s) const { return new Mul(a->substitute(s), b->substitute(s)); }
    Expression* derivative(const Variable& var) const;
    Expression* simplify() const;
    Mul* construct(Expression* _a, Expression* _b) const { return new Mul(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") * (" + b->toString() + ")"; }
};
struct Div : public BinaryOp {
    Div(Expression* _a, Expression* _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return a->evaluate() / b->evaluate(); }
    Number evaluate(Number _a, Number _b) const { return _a / _b; }
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& c) const { return c.CreateFDiv(a->evaluate(c), b->evaluate(c)); }
#endif
    Expression* substitute(const Subst& s) const { return new Div(a->substitute(s), b->substitute(s)); }
    Expression* derivative(const Variable& var) const;
    Expression* simplify() const;
    Div* construct(Expression* _a, Expression* _b) const { return new Div(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") / (" + b->toString() + ")"; }
};
struct Pow : public BinaryOp {
    Pow(Expression* _a, Expression* _b) : BinaryOp(_a, _b) { }
    Number evaluate() const { return std::pow(a->evaluate(), b->evaluate()); }
    Number evaluate(Number _a, Number _b) const { return std::pow(_a, _b); }
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& c) const {
        AsmJit::XMMVar ret(c.newXMM());
        Number (*pow)(Number, Number) = &std::pow;
        AsmJit::ECall* ctx = c.call((void*)pow);
        ctx->setPrototype(AsmJit::CALL_CONV_DEFAULT, AsmJit::FunctionBuilder2<Number, Number, Number>());
        ctx->setArgument(0, a->evaluate(c));
        ctx->setArgument(1, b->evaluate(c));
        ctx->setReturn(ret);
        return ret;
    }
#endif
    Expression* substitute(const Subst& s) const { return new Pow(a->substitute(s), b->substitute(s)); }
    Expression* derivative(const Variable& var) const;
    Expression* simplify() const;
    Pow* construct(Expression* _a, Expression* _b) const { return new Pow(_a, _b); }
    std::string toString() const { return std::string("(") + a->toString() + ") ^ (" + b->toString() + ")"; }
};
struct PowInt : public UnaryOp {
    int b;
    PowInt(Expression* _a, int _b) : UnaryOp(_a), b(_b) { }
    Number evaluate() const;
    Number evaluate(Number _a) const;
    Vector evaluateVector(std::size_t size) const;
#ifdef USE_LLVM
    llvm::Value* evaluate(llvm::IRBuilder<>& c) const {
        llvm::Value* _a = a->evaluate(c);
        llvm::Value* ret = llvm::ConstantFP::get(llvm::getGlobalContext(), 1.);
        int _b = b;
        bool neg = b < 0;
        if (neg) _b = -_b;
        while (_b > 0) {
            if (_b & 1) ret = c.CreateFMul(ret, _a);
            _b >>= 1;
            _a = c.CreateFMul(_a, _a);
        }
        if (neg) {
            return c.CreateFDiv(llvm::ConstantFP::get(llvm::getGlobalContext(), 1.), ret);
        }
        return ret;
    }
#endif
    Expression* substitute(const Subst& s) const { return new PowInt(a->substitute(s), b); }
    Expression* derivative(const Variable& var) const;
    Expression* simplify() const;
    PowInt* construct(Expression* _a) const { return new PowInt(_a, b); }
    std::string toString() const { return std::string("(") + a->toString() + ") ^ (!" + boost::lexical_cast<std::string>(b) + ")"; }
};
struct PolyGamma : public UnaryOp {
    int b;
    PolyGamma(Expression* _a, int _b) : UnaryOp(_a), b(_b) { }
    Number evaluate() const;
    Number evaluate(Number _a) const;
    Vector evaluateVector(std::size_t size) const;
    Expression* substitute(const Subst& s) const { return new PolyGamma(a->substitute(s), b); }
    Expression* derivative(const Variable& var) const;
    PolyGamma* construct(Expression* _a) const { return new PolyGamma(_a, b); }
    std::string toString() const { return std::string("psi_") + boost::lexical_cast<std::string>(b) + "(" + a->toString() + ")"; }
};
struct Equation : public Thing {
    Expression* a, * b;
    Equation(Expression* _a, Expression* _b) : a(_a), b(_b) { }
    ~Equation() { delete a; delete b; }
    std::string toString() const { return a->toString() + " = " + b->toString(); }
};
struct Function : public Thing {
	struct ArgumentMismatchException { };
	struct UnboundVariableException { };
	std::vector<Variable> bound_args;
	Expression* body;
	Function() : body(NULL) { }
	Function(std::vector<Variable> _args, Expression* _body) : bound_args(_args), body(_body) { check(); }
	~Function() { delete body; }
	void check() {
		std::set<Variable> vars;
		body->variables(vars);
		if (vars.size() != bound_args.size()) throw UnboundVariableException();
        for (std::size_t i = 0; i < bound_args.size(); ++i) {
            if (vars.find(bound_args[i]) == vars.end()) throw UnboundVariableException();
        }
	}
	Expression* call(const std::vector<Expression*>& args) {
		if (args.size() != bound_args.size()) throw ArgumentMismatchException();
		Expression::Subst s;
		for (std::size_t i = 0; i < bound_args.size(); ++i) {
			s.insert(std::make_pair(bound_args[i], args[i]));
		}
		return body->substitute(s);
	}
};

#endif
