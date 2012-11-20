cmake_minimum_required(VERSION 2.8)

project(grapher)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(Qt4 REQUIRED)
find_package(GSL REQUIRED)
find_package(Vc REQUIRED)

set(asmjit_SOURCES
    AsmJit/AssemblerX86X64.cpp
    AsmJit/CodeGenerator.cpp
    AsmJit/Compiler.cpp
    AsmJit/CompilerX86X64.cpp
    AsmJit/CpuInfo.cpp
    AsmJit/Defs.cpp
    AsmJit/DefsX86X64.cpp
    AsmJit/Logger.cpp
    AsmJit/MemoryManager.cpp
    AsmJit/MemoryMarker.cpp
    AsmJit/OperandX86X64.cpp
    AsmJit/Platform.cpp
    AsmJit/Util.cpp)

set(grapher_SOURCES
    main.cpp
    sse_mathfun.cpp
    matrix_inverse.cpp
    Expression.cpp
    Expression_jit.cpp
    Graph.cpp
    Graph2D.cpp
    Graph3D.cpp
    Grapher2D.cpp
    Grapher3D.cpp
    Parser.cpp
    MainWindow.cpp
    Window2DSettings.cpp
    Window3DSettings.cpp
    Graph2DProperties.cpp
    Graph3DProperties.cpp
    UI2DFactory.cpp
    UI3DFactory.cpp
    UI.cpp
    Render3D.cpp
    ${asmjit_SOURCES}
)
set(grapher_HEADERS
    Expression.h
    Graph.h
    Graph2D.h
    Graph3D.h
    Grapher2D.h
    Grapher3D.h
    Parser.h
    MainWindow.h
    Window2DSettings.h
    Window3DSettings.h
    Graph2DProperties.h
    Graph3DProperties.h
    UIFactory.h
    UI2DFactory.h
    UI3DFactory.h
    UI.h
    Render3D.h
)
set(grapher_UI
    MainWindow.ui
    UI.ui
    Window2DSettings.ui
    Window3DSettings.ui
    Graph2DProperties.ui
    Graph3DProperties.ui
)

qt4_wrap_cpp(grapher_HEADERS_MOC ${grapher_HEADERS})
qt4_wrap_ui(grapher_SOURCES_UI ${grapher_UI})

include_directories(${GSL_INCLUDE_DIRS})
include_directories(.)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${Vc_INCLUDE_DIR})

include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})
add_definitions(${Vc_DEFINITIONS})
add_definitions(-std=gnu++0x)

add_executable(grapher ${grapher_SOURCES} ${grapher_HEADERS_MOC} ${grapher_SOURCES_UI})
target_link_libraries(grapher ${QT_LIBRARIES} ${GSL_LIBRARIES} ${Vc_LIBRARIES})
