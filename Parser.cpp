#include "Parser.h"

#include "Expression.h"
#include <list>
#include <stack>
#include <iostream>
#include <sstream>
#include <boost/variant.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/unordered_map.hpp>

bool isIntegral(Number a) {
    const Number b = floor(a + 0.5);
    const Number epsilon = 1e-4;
    return (a < b + epsilon) && (b < a + epsilon);
}

int rnd(Number a) {
    return static_cast<int>(a + 0.5);
}

struct printer : public boost::static_visitor<> {
    void operator()(std::string s) const {
        std::cerr << s << std::endl;
    }
    void operator()(char s) const {
        std::cerr << s << std::endl;
    }
    void operator()(Number s) const {
        std::cerr << s << std::endl;
    }
};
struct shunt_visitor : public boost::static_visitor<> {
    typedef boost::variant<std::string, char> optoken;
    std::stack<Thing*> output;
    std::stack<optoken> operators;
    std::stack<int> argcounts;
    const boost::unordered_map<std::string, Expression*>& variables;
    shunt_visitor(const boost::unordered_map<std::string, Expression*>& _v) : variables(_v) { }
    ~shunt_visitor() {
        while (!output.empty()) {
            delete output.top();
            output.pop();
        }
    }
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
    void operator()(std::string name) {
        boost::unordered_map<std::string, Expression*>::const_iterator it = variables.find(name);
        if (it != variables.end()) {
            output.push(it->second->copy());
        } else {
            // a function?
            operators.push(name);
            argcounts.push(1);
        }
    }
    Expression* pop_expression() {
        if (output.empty()) throw Parser::UnexpectedOperatorException();
        Thing* top = output.top();
        output.pop();
        Expression* ret = dynamic_cast<Expression*>(top);
        if (ret == NULL) throw Parser::ExpressionTypeException();
        return ret;
    }
    void apply_operator(optoken op) {
        if (std::string* _func = boost::get<std::string>(&op)) {
            std::string& func = *_func;
            if (argcounts.empty()) throw Parser::ArgumentCountException("before " + func);
            int argc = argcounts.top();
            argcounts.pop();
#define UNARY_FUNCTION(class, name) \
if (func == #name) { \
    if (argc != 1) throw Parser::ArgumentCountException(func); \
    output.push(new class(pop_expression())); \
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
                Expression* z = pop_expression();
                Expression* n = pop_expression();
                Constant* c = dynamic_cast<Constant*>(n);
                if (c == NULL || !isIntegral(c->c)) throw Parser::ExpressionTypeException();
                output.push(new PolyGamma(z, rnd(c->c)));
            } else {
                throw Parser::UnknownFunctionException(func);
            }
#undef UNARY_FUNCTION
        } else if (char* _c = boost::get<char>(&op)) {
            char c = *_c;
            Expression* b = pop_expression();
            Expression* a = pop_expression();
            switch (c) {
                case '*': output.push(new Mul(a, b)); break;
                case '/': output.push(new Div(a, b)); break;
                case '+': output.push(new Add(a, b)); break;
                case '-': output.push(new Sub(a, b)); break;
                case '^': {
                    Constant* bc = dynamic_cast<Constant*>(b);
                    if (bc != NULL && isIntegral(bc->c)) {
                        output.push(new PowInt(a, rnd(bc->c)));
                        delete bc;
                    } else {
                        output.push(new Pow(a, b));
                    }
                    break;
                }
                case '=': output.push(new Equation(a, b)); break;
            }
        } else {
            throw "EVERYTHING'S BAD";
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
        output.push(new Constant(num));
    }
    Thing* finish() {
        while (!operators.empty()) pop_apply();
        if (output.size() != 1) throw Parser::MissingOperatorException();
        Thing* ret = output.top();
        output.pop();
        return ret;
    }
};

Thing* Parser::parse(const std::string& str, const boost::unordered_map<std::string, Expression*>& variables) {
    assert(isIntegral(1.0));
    assert(isIntegral(2.0));
    assert(!isIntegral(2.1));
    assert(!isIntegral(515.02));
    std::cerr << str << std::endl;
    typedef boost::variant<std::string, char, Number> token;
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
            tokens.push_back(name);
        } else {
            tokens.push_back(c);
        }
    }
    for (std::list<token>::iterator next = tokens.begin(), p = next++; next != tokens.end(); p = next++) {
        if (boost::get<char>(&*p) == NULL && boost::get<char>(&*next) == NULL) {
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