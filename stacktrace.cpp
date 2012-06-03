#include "stacktrace.h"

#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#include <dbghelp.h>

void printStack() {
    unsigned int i;
    void* stack[100];
    unsigned short frames;
    SYMBOL_INFO* symbol;
    HANDLE process;

    process = GetCurrentProcess();

    SymInitialize(process, NULL, TRUE);

    frames = CaptureStackBackTrace(0, 100, stack, NULL);
    symbol = (SYMBOL_INFO *)calloc(sizeof(SYMBOL_INFO) + 256 * sizeof(char), 1);
    symbol->MaxNameLen = 255;
    symbol->SizeOfStruct = sizeof(SYMBOL_INFO);

    printf("%d frames\n", (int)frames);

    for (i = 0; i < frames; i++) {
        SymFromAddr(process, (DWORD64)(stack[i]), 0, symbol);

        printf("%i: %s - 0x%0X\n", frames - i - 1, symbol->Name, symbol->Address);
    }

    free(symbol);
    fflush(stdout);
}
#else
#include <execinfo.h>
#include <cxxabi.h>


void printStack() {
    printf("TODO: stack trace");
}

#endif
