cmake_minimum_required(VERSION 2.8)

project(grapher)

find_package(Qt4 REQUIRED)
find_package(GSL REQUIRED)

set(grapher_SOURCES main.cpp Expression.cpp Graph.cpp Grapher.cpp Parser.cpp MainWindow.cpp WindowSettings.cpp GraphProperties.cpp)
set(grapher_HEADERS Expression.h Graph.h Grapher.h Parser.h MainWindow.h WindowSettings.h GraphProperties.h)
set(grapher_UI MainWindow.ui WindowSettings.ui GraphProperties.ui)

# set(asmjit_SOURCES AsmJit/Assembler.cpp AsmJit/AssemblerX86X64.cpp AsmJit/CodeGenerator.cpp AsmJit/Compiler.cpp AsmJit/CompilerX86X64.cpp AsmJit/CpuInfo.cpp AsmJit/Defs.cpp AsmJit/DefsX86X64.cpp AsmJit/Logger.cpp AsmJit/MemoryManager.cpp AsmJit/Operand.cpp AsmJit/OperandX86X64.cpp AsmJit/Platform.cpp AsmJit/Util.cpp)
# set(asmjit_HEADERS AsmJit/ApiBegin.h AsmJit/ApiEnd.h AsmJit/AsmJit.h AsmJit/Assembler.h AsmJit/AssemblerX86X64.h AsmJit/Build.h AsmJit/CodeGenerator.h AsmJit/Compiler.h AsmJit/CompilerX86X64.h AsmJit/Config.h AsmJit/CpuInfo.h AsmJit/Defs.h AsmJit/DefsX86X64.h AsmJit/Logger.h AsmJit/MemoryManager.h AsmJit/Operand.h AsmJit/OperandX86X64.h AsmJit/Platform.h AsmJit/Util.h AsmJit/Util_p.h)
# include_directories(AsmJit)

qt4_wrap_cpp(grapher_HEADERS_MOC ${grapher_HEADERS})
qt4_wrap_ui(grapher_SOURCES_UI ${grapher_UI})

include_directories(${GSL_INCLUDE_DIRS})
include_directories(.)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

add_definitions(-std=c++0x)

add_executable(grapher ${grapher_SOURCES} ${grapher_HEADERS_MOC} ${grapher_SOURCES_UI})
target_link_libraries(grapher ${QT_LIBRARIES} ${GSL_LIBRARIES})
