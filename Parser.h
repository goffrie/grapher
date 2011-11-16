#ifndef _PARSER_H_
#define _PARSER_H_

#include <string>
#include <exception>
#include <unordered_map>
#include <memory>

struct Thing;
struct Expression;

struct Parser {
    struct InvalidNumberException : public std::exception {
        const char* what() const throw() {
            return "Unparseable number";
        }
    };
    struct UnknownOperatorException : public std::exception {
        std::string msg;
        UnknownOperatorException(char op) : msg("Unknown operator: ") { msg += op; }
        ~UnknownOperatorException() throw() { }
        const char* what() const throw() {
            return msg.c_str();
        }
    };
    struct UnknownFunctionException : public std::exception {
        std::string msg;
        UnknownFunctionException(std::string func) : msg("Unknown function: ") { msg += func; }
        ~UnknownFunctionException() throw() { }
        const char* what() const throw() {
            return msg.c_str();
        }
    };
    struct InvalidCommaException : public std::exception {
        const char* what() const throw() {
            return "Invalid comma!";
        }
    };
    struct MismatchedParenthesesException : public std::exception {
        const char* what() const throw() {
            return "Mismatched parentheses!";
        }
    };
    struct ArgumentCountException : public std::exception {
        std::string msg;
        ArgumentCountException(std::string func) : msg("Wrong number of arguments to function: ") { msg += func; }
        ~ArgumentCountException() throw() { }
        const char* what() const throw() {
            return msg.c_str();
        }
    };
    struct ExpressionTypeException : public std::exception {
        const char* what() const throw() {
            return "Wrong type of expression!";
        }
    };
    struct UnexpectedOperatorException : public std::exception {
        const char* what() const throw() {
            return "Didn't expect that operator here!";
        }
    };
    struct MissingOperatorException : public std::exception {
        const char* what() const throw() {
            return "Expected more operators than that!";
        }
    };
    static std::unique_ptr<Thing> parse(const std::string& str, const std::unordered_map<std::string, Expression*>& variables);
};


#endif
