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
        std::cerr << s << ' ';
    }
    template<class T> void operator()(const Name<T>& s) const {
        std::cerr << (std::string&)s << ' ';
    }
    void operator()(char s) const {
        std::cerr << s << ' ';
    }
    void operator()(Number s) const {
        std::cerr << s << ' ';
    }
};
struct shunt_visitor : public boost::static_visitor<> {
    typedef boost::variant<Name<FunctionTag>, char> optoken;
    std::vector<std::unique_ptr<Thing> > output;
    std::vector<optoken> operators;
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
            case '^': case '_': return RightAssociative;
            default: throw Parser::UnknownOperatorException(op);
        }
    }
    static int precedence(char op) {
        switch (op) {
            case '^': return 10;
            case '*': case '/': return 6;
            case '_': return 4;
            case '+': case '-': return 3;
            case '=': return 0;
            default: throw Parser::UnknownOperatorException(op);
        }
    }
    static bool unary(char op) {
        return op == '_';
    }
    void operator()(const Name<FunctionTag>& name) {
        operators.push_back(name);
        argcounts.push(1);
    }
    void operator()(const Name<VariableTag>& name) {
        output.push_back(variables.find(name)->second->ecopy());
    }
    EPtr pop_expression() {
        if (output.empty()) throw Parser::UnexpectedOperatorException();
        // std::cerr << output.back()->toString() << std::endl;
        EPtr ret(dynamic_unique_cast<Expression>(std::move(output.back())));
        output.pop_back();
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
    output.push_back(class::create(pop_expression())); \
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
                output.push_back(PolyGamma::create(std::move(z), rnd(c->c)));
            } else {
                throw Parser::UnknownFunctionException(func);
            }
#undef UNARY_FUNCTION
        } else if (char* _c = boost::get<char>(&op)) {
            char c = *_c;
            if (c == '_') {
                output.push_back(Neg::create(pop_expression()));
            } else {
                EPtr b(pop_expression());
                EPtr a(pop_expression());
                switch (c) {
                    case '*': output.push_back(Mul::create(std::move(a), std::move(b))); break;
                    case '/': output.push_back(Div::create(std::move(a), std::move(b))); break;
                    case '+': output.push_back(Add::create(std::move(a), std::move(b))); break;
                    case '-': output.push_back(Sub::create(std::move(a), std::move(b))); break;
                    case '^': {
                        Constant* bc = dynamic_cast<Constant*>(b.get());
                        if (bc != NULL && isIntegral(bc->c)) {
                            output.push_back(PowInt::create(std::move(a), rnd(bc->c)));
                        } else {
                            output.push_back(Pow::create(std::move(a), std::move(b)));
                        }
                        break;
                    }
                    case '=': output.push_back(Equation::create(std::move(a), std::move(b))); break;
                    default: std::terminate();
                }
            }
        } else {
            qDebug() << "EVERYTHING'S BAD";
            std::terminate();
        }
    }
    void pop_apply() {
        apply_operator(operators.back());
        operators.pop_back();
    }
    void operator()(char op) {
        // std::cerr << op << std::endl;
        switch (op) {
        case ',': {
            char *c;
            while (((c = NULL), (!operators.empty())) && (c = boost::get<char>(&operators.back())) != NULL && *c != '(') {
                pop_apply();
            }
            if (c == NULL) throw Parser::InvalidCommaException();
            int n = argcounts.top();
            argcounts.pop();
            argcounts.push(n+1);
            break;
        }
        case '_':
            operators.push_back(op);
            break;
        case '*': case '/': case '-': case '+': case '^': case '=': {
            int last = operators.size();
            for (int i = last; i--; ) {
                char *_c = boost::get<char>(&operators[i]);
                if (_c == NULL) throw Parser::UnexpectedOperatorException();
                char c = *_c;
                if (c == '(') break;
                if (
                    (associativity(op) == LeftAssociative) ?
                        (precedence(op) <= precedence(c)) :
                        (precedence(op) <  precedence(c))
                ) {
                    last = i;
                }
            }
            while (operators.size() > last) pop_apply();
            operators.push_back(op);
            break;
        }
        case '(':
            operators.push_back(op);
            break;
        case ')':
            char *c;
            while (((c = NULL), (!operators.empty())) && (c = boost::get<char>(&operators.back())) != NULL && *c != '(') {
                pop_apply();
            }
            if (c == NULL) throw Parser::MismatchedParenthesesException();
            assert(boost::get<char>(operators.back()) == '(');
            operators.pop_back();
            if (!operators.empty() && boost::get<Name<FunctionTag>>(&operators.back()) != NULL) {
                pop_apply(); // apply function
            }
            break;
        default:
            throw Parser::UnknownOperatorException(op);
        }
    }
    void operator()(Number num) {
        output.push_back(Constant::create(num));
    }
    std::unique_ptr<Thing> finish() {
        while (!operators.empty()) pop_apply();
        if (output.size() != 1) throw Parser::MissingOperatorException();
        std::unique_ptr<Thing> ret(std::move(output.back()));
        output.pop_back();
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
    for (std::list<token>::iterator next = tokens.begin(), p; next != tokens.end(); p = next++) {
        char* a = (next == tokens.begin()) ? NULL : boost::get<char>(&*p),
            * b = boost::get<char>(&*next);
        if (next != tokens.begin() && (!a || *a == ')') && (!b || (*b == '(' && !boost::get<Name<FunctionTag>>(&*p)))) {
            // implicit multiplication
            next = tokens.insert(next, '*');
        } else if ((next == tokens.begin() || (a && *a != ')')) && b && *b == '-') {
            // unary minus
            *b = '_';
        }
    }
    shunt_visitor v(variables);
/*    for (std::list<token>::iterator p = tokens.begin(); p != tokens.end(); ++p) {
        boost::apply_visitor(printer(), *p);
    }
    std::cerr << std::endl;*/
    for (std::list<token>::iterator p = tokens.begin(); p != tokens.end(); ++p) {
        boost::apply_visitor(v, *p);
/*        for (int i = 0; i < v.output.size(); ++i) {
            std::cerr << v.output[i]->toString() << " | ";
        }
        for (int i = v.operators.size(); i--; ) {
            boost::apply_visitor(printer(), v.operators[i]);
        }
        std::cerr << std::endl;*/
    }
    return v.finish();
}