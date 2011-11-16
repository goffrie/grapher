#include "Parser.h"

#include "Expression.h"
#include <list>
#include <stack>
#include <iostream>
#include <sstream>
#include <memory>
#include <boost/variant.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/unordered_map.hpp>
#include <QDebug>

#include "dynamic_unique_cast.h"
#include "util.h"

struct FunctionTag {
};
struct VariableTag {
};
template<class Tag> struct Name : public std::string {
    explicit Name(const std::string& s) : std::string(s) { }
};
/*
struct VariableName : public std::string {
    explicit VariableName(const std::string& s) : std::string(s) { }
};

template<class T> struct FunctionBuilder {
};
*/

struct printer : public boost::static_visitor<> {
    void operator()(const std::string& s) const {
        std::cerr << s << std::endl;
    }
    template<class T> void operator()(const Name<T>& s) const {
        std::cerr << (std::string&)s << std::endl;
    }
    void operator()(char s) const {
        std::cerr << s << std::endl;
    }
    void operator()(Number s) const {
        std::cerr << s << std::endl;
    }
};
struct shunt_visitor : public boost::static_visitor<> {
    typedef boost::variant<Name<FunctionTag>, char> optoken;
    std::stack<std::unique_ptr<Thing> > output;
    std::stack<optoken> operators;
    std::stack<int> argcounts;
    const boost::unordered_map<std::string, Expression*>& variables;
    shunt_visitor(const boost::unordered_map<std::string, Expression*>& _v) : variables(_v) { }
    enum Associativity {
        LeftAssociative,
        RightAssociative
    };
    static Associativity associativity(char op) {
        switch (op) {
            case '*': case '/': case '-': case '+': case '=': return LeftAssociative;
            case '^': return RightAssociative;
            default: throw Parser::UnknownOperatorException(op);
        }
    }
    static int precedence(char op) {
        switch (op) {
            case '^': return 3;
            case '*': case '/': return 2;
            case '+': case '-': return 1;
            case '=': return 0;
            default: throw Parser::UnknownOperatorException(op);
        }
    }
    void operator()(const Name<FunctionTag>& name) {
        operators.push(name);
        argcounts.push(1);
    }
    void operator()(const Name<VariableTag>& name) {
        output.push(variables.find(name)->second->ecopy());
    }
    EPtr pop_expression() {
        if (output.empty()) throw Parser::UnexpectedOperatorException();
        EPtr ret(dynamic_unique_cast<Expression>(std::move(output.top())));
        output.pop();
        if (!ret) throw Parser::ExpressionTypeException();
        return std::move(ret);
    }
    void apply_operator(optoken op) {
        if (Name<FunctionTag>* _func = boost::get<Name<FunctionTag> >(&op)) {
            std::string& func = *_func;
            if (argcounts.empty()) throw Parser::ArgumentCountException("before " + func);
            int argc = argcounts.top();
            argcounts.pop();
#define UNARY_FUNCTION(class, name) \
if (func == #name) { \
    if (argc != 1) throw Parser::ArgumentCountException(func); \
    output.push(class::create(pop_expression())); \
} else
            UNARY_FUNCTION(Exp, exp)
            UNARY_FUNCTION(Ln, ln)
            UNARY_FUNCTION(Sqrt, sqrt)
            UNARY_FUNCTION(Sin, sin)
            UNARY_FUNCTION(Cos, cos)
            UNARY_FUNCTION(Tan, tan)
            UNARY_FUNCTION(Asin, asin)
            UNARY_FUNCTION(Acos, acos)
            UNARY_FUNCTION(Atan, atan)
            UNARY_FUNCTION(Gamma, gamma)
            /* else */ if (func == "psi") {
                if (argc != 2) throw Parser::ArgumentCountException(func);
                EPtr z(pop_expression());
                std::unique_ptr<Constant> c = dynamic_unique_cast<Constant>(pop_expression());
                if (!c || !isIntegral(c->c)) throw Parser::ExpressionTypeException();
                output.push(PolyGamma::create(std::move(z), rnd(c->c)));
            } else {
                throw Parser::UnknownFunctionException(func);
            }
#undef UNARY_FUNCTION
        } else if (char* _c = boost::get<char>(&op)) {
            char c = *_c;
            EPtr b(pop_expression());
            EPtr a(pop_expression());
            switch (c) {
                case '*': output.push(Mul::create(std::move(a), std::move(b))); break;
                case '/': output.push(Div::create(std::move(a), std::move(b))); break;
                case '+': output.push(Add::create(std::move(a), std::move(b))); break;
                case '-': output.push(Sub::create(std::move(a), std::move(b))); break;
                case '^': {
                    Constant* bc = dynamic_cast<Constant*>(b.get());
                    if (bc != NULL && isIntegral(bc->c)) {
                        output.push(PowInt::create(std::move(a), rnd(bc->c)));
                    } else {
                        output.push(Pow::create(std::move(a), std::move(b)));
                    }
                    break;
                }
                case '=': output.push(Equation::create(std::move(a), std::move(b))); break;
                default: std::terminate();
            }
        } else {
            qDebug() << "EVERYTHING'S BAD";
            std::terminate();
        }
    }
    void pop_apply() {
        apply_operator(operators.top());
        operators.pop();
    }
    void operator()(char op) {
        std::cerr << op << std::endl;
        switch (op) {
        case ',': {
            char *c;
            while (((c = NULL), (!operators.empty())) && (c = boost::get<char>(&operators.top())) != NULL && *c != '(') {
                pop_apply();
            }
            if (c == NULL) throw Parser::InvalidCommaException();
            int n = argcounts.top();
            argcounts.pop();
            argcounts.push(n+1);
            break;
        }
        case '*': case '/': case '-': case '+': case '^': case '=': {
            char* c;
            while (!operators.empty() && (c = boost::get<char>(&operators.top())) != NULL &&
                (*c != '(') &&
                ((associativity(op) == LeftAssociative) ?
                    (precedence(op) <= precedence(*c)) :
                    (precedence(op) < precedence(*c)))
            ) {
                pop_apply();
            }
            operators.push(op);
            break;
        }
        case '(':
            operators.push(op);
            break;
        case ')':
            char *c;
            while (((c = NULL), (!operators.empty())) && (c = boost::get<char>(&operators.top())) != NULL && *c != '(') {
                pop_apply();
            }
            if (c == NULL) throw Parser::MismatchedParenthesesException();
            assert(boost::get<char>(operators.top()) == '(');
            operators.pop();
            if (!operators.empty() && boost::get<std::string>(&operators.top()) != NULL) {
                pop_apply(); // apply function
            }
            break;
        default:
            throw Parser::UnknownOperatorException(op);
        }
    }
    void operator()(Number num) {
        output.push(Constant::create(num));
    }
    std::unique_ptr<Thing> finish() {
        while (!operators.empty()) pop_apply();
        if (output.size() != 1) throw Parser::MissingOperatorException();
        std::unique_ptr<Thing> ret(std::move(output.top()));
        output.pop();
        return std::move(ret);
    }
};

std::unique_ptr<Thing> Parser::parse(const std::string& str, const boost::unordered_map<std::string, Expression*>& variables) {
    typedef boost::variant<Name<FunctionTag>, Name<VariableTag>, char, Number> token;
    std::list<token> tokens;
    for (std::size_t i = 0; i < str.size(); ++i) {
        char c = str[i];
        if (std::isspace(c)) continue;
        if (std::isdigit(c) || c == '.') {
            std::string num(1, c);
            while ((i+1) < str.size() && (std::isdigit(str[i+1]) || str[i+1] == '.')) {
                num += str[++i];
            }
            std::istringstream stream(num);
            Number n;
            if (!(stream >> n)) throw InvalidNumberException();
            tokens.push_back(n);
        } else if (std::isalpha(c)) {
            std::string name(1, c);
            while ((i+1) < str.size() && std::isalnum(str[i+1])) {
                name += str[++i];
            }
            if (variables.find(name) != variables.end()) {
                tokens.push_back(Name<VariableTag>(name));
            } else {
                tokens.push_back(Name<FunctionTag>(name));
            }
        } else {
            tokens.push_back(c);
        }
    }
    for (std::list<token>::iterator next = tokens.begin(), p = next++; next != tokens.end(); p = next++) {
        char* a = boost::get<char>(&*p),
            * b = boost::get<char>(&*next);
        if ((!a || *a == ')') && (!b || (*b == '(' && !boost::get<Name<FunctionTag>>(&*p)))) {
            next = tokens.insert(next, '*');
        }
    }
    shunt_visitor v(variables);
    for (std::list<token>::iterator p = tokens.begin(); p != tokens.end(); ++p) {
        boost::apply_visitor(printer(), *p);
        boost::apply_visitor(v, *p);
    }
    return v.finish();
}