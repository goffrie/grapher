#include "Expression.h"

// TODO: make this not terrible
#undef V

#include <llvm/IR/DerivedTypes.h>
#include <llvm/ExecutionEngine/ExecutionEngine.h>
#include <llvm/ExecutionEngine/JIT.h>
#include <llvm/IR/IRBuilder.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/IR/Module.h>
#include <llvm/Analysis/Verifier.h>
#include <llvm/Analysis/Passes.h>
#include <llvm/IR/DataLayout.h>
#include <llvm/Transforms/Scalar.h>
#include <llvm/Support/TargetSelect.h>
#include <llvm/IR/Intrinsics.h>

using namespace llvm;

void dispose(Function* f) {
    if (f) f->eraseFromParent();
}

MathContext MathContext::defaultContext() {
    static bool herp = true;
    if (herp) { herp = false; InitializeNativeTarget(); }
    Module* module = new Module("derp", getGlobalContext());
    ExecutionEngine* engine = EngineBuilder(module).create();
    
    FunctionPassManager* fpm = new FunctionPassManager(module);
    fpm->add(new DataLayout(*engine->getDataLayout()));
    fpm->add(createBasicAliasAnalysisPass());
    fpm->add(createInstructionCombiningPass());
    fpm->add(createReassociatePass());
    fpm->add(createGVNPass());
    fpm->add(createCFGSimplificationPass());
    fpm->doInitialization();

    return { engine, fpm, module };
}

WEvalFunc Expression::evaluator(MathContext& context) const {
    FunctionType* type = FunctionType::get(Type::getDoubleTy(getGlobalContext()), false);
    Function* func = Function::Create(type, Function::ExternalLinkage, "", context.module);
    BasicBlock* block = BasicBlock::Create(getGlobalContext(), "entrypoint", func);
    auto builder = IRBuilder<>(block);
    builder.CreateRet(evaluateJIT(builder, context));
    verifyFunction(*func);
    func->dump();
    context.fpm->run(*func); // optimize
    func->dump();
    // Run JIT
    auto fptr = reinterpret_cast<EvalFunc>(context.jitEngine->getPointerToFunction(func));
    return WEvalFunc(func, fptr);
}

Value* Variable::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    switch (id->type) {
        case Id::Constant:
        case Id::Vector:
            return builder.CreateCast(Instruction::FPExt, 
                builder.CreateLoad(
                    llvm::Constant::getIntegerValue(Type::getFloatPtrTy(builder.getContext()),
                    APInt(sizeof(float*)*8, reinterpret_cast<uint64_t>(id->p)))),
                Type::getDoubleTy(builder.getContext()));
        default:
            throw this;
    }
}

Value* ::Constant::evaluateJIT(IRBuilder<>& builder, MathContext&) const {
    return ConstantFP::get(Type::getDoubleTy(builder.getContext()), c);
}

Value* Neg::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    return builder.CreateFNeg(a->evaluateJIT(builder, context));
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

#define UNARY_INTRINSIC_FUNCTION(name, func) \
Value* name::evaluateJIT(IRBuilder<>& builder, MathContext& context) const { \
    return builder.CreateCall(Intrinsic::getDeclaration(context.module, Intrinsic::func, Type::getDoubleTy(builder.getContext())), \
        a->evaluateJIT(builder, context)); \
}

UNARY_INTRINSIC_FUNCTION(Exp, exp)
UNARY_INTRINSIC_FUNCTION(Ln, log)
UNARY_INTRINSIC_FUNCTION(Sqrt, sqrt)
UNARY_INTRINSIC_FUNCTION(Sin, sin)
UNARY_INTRINSIC_FUNCTION(Cos, cos)

Value* Tan::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    Value* theta = a->evaluateJIT(builder, context);
    return builder.CreateFDiv(
        builder.CreateCall(Intrinsic::getDeclaration(context.module, Intrinsic::sin, Type::getDoubleTy(builder.getContext())), theta),
        builder.CreateCall(Intrinsic::getDeclaration(context.module, Intrinsic::cos, Type::getDoubleTy(builder.getContext())), theta));
}

#define UNARY_EXTERNAL_FUNCTION(name, func) \
Value* name::evaluateJIT(IRBuilder<>& builder, MathContext& context) const { \
    return builder.CreateCall(getUnaryGlobalFunction(builder, context, #func), a->evaluateJIT(builder, context)); \
}

static Function* getUnaryGlobalFunction(IRBuilder<>& builder, MathContext& context, const char* name) {
    Function* exists = context.module->getFunction(name);
    if (exists) return exists;
    std::vector<Type*> argtypes(1, Type::getDoubleTy(builder.getContext()));
    FunctionType* type = FunctionType::get(Type::getDoubleTy(builder.getContext()), argtypes, false);
    return Function::Create(type, Function::ExternalLinkage, name, context.module);
}

UNARY_EXTERNAL_FUNCTION(Asin, asin)
UNARY_EXTERNAL_FUNCTION(Acos, acos)
UNARY_EXTERNAL_FUNCTION(Atan, atan)
UNARY_EXTERNAL_FUNCTION(Gamma, gsl_sf_gamma)

Value* Add::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    return builder.CreateFAdd(a->evaluateJIT(builder, context), b->evaluateJIT(builder, context));
}

Value* Sub::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    return builder.CreateFSub(a->evaluateJIT(builder, context), b->evaluateJIT(builder, context));
}

Value* Mul::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    return builder.CreateFMul(a->evaluateJIT(builder, context), b->evaluateJIT(builder, context));
}

Value* Div::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    return builder.CreateFDiv(a->evaluateJIT(builder, context), b->evaluateJIT(builder, context));
}

Value* PowInt::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    return builder.CreateCall2(Intrinsic::getDeclaration(context.module, Intrinsic::powi, Type::getDoubleTy(builder.getContext())),
                               a->evaluateJIT(builder, context),
                               ConstantInt::getSigned(Type::getInt32Ty(builder.getContext()), b));
}

Value* Polynomial::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    if (!left.get()) return right->evaluateJIT(builder, context);
    return builder.CreateFAdd(
        builder.CreateFMul(
            left->evaluateJIT(builder, context),
            var.evaluateJIT(builder, context)
        ),
        right->evaluateJIT(builder, context)
    );
}

Value* Pow::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    return builder.CreateCall2(Intrinsic::getDeclaration(context.module, Intrinsic::pow, Type::getDoubleTy(builder.getContext())),
                               a->evaluateJIT(builder, context),
                               b->evaluateJIT(builder, context));
}

Value* PolyGamma::evaluateJIT(IRBuilder<>& builder, MathContext& context) const {
    Function* func = context.module->getFunction("gsl_sf_psi_n");
    if (!func) {
        std::vector<Type*> argtypes;
        argtypes.push_back(Type::getInt32Ty(builder.getContext()));
        argtypes.push_back(Type::getDoubleTy(builder.getContext()));
        FunctionType* type = FunctionType::get(Type::getDoubleTy(builder.getContext()), argtypes, false);
        func = Function::Create(type, Function::ExternalLinkage, "gsl_sf_psi_n", context.module);
    }
    return builder.CreateCall2(func, ConstantInt::getSigned(Type::getInt32Ty(builder.getContext()), b),
        a->evaluateJIT(builder, context)
    );
}
