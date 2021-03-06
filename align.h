#ifndef ALIGN_H
#define ALIGN_H

#include <cstdlib>
#include <cstdint>
#include <new>
#include <Vc/Vc>

#define ALIGN_AMT (VC_FLOAT_V_SIZE*sizeof(float))

inline void* aligned_malloc(size_t size) {
    void* pa = std::malloc(size + (ALIGN_AMT-1) + sizeof(void *));
    if (!pa) throw std::bad_alloc();
    void* ptr = (void*) ( ((uintptr_t)pa + sizeof(void *) + (ALIGN_AMT-1)) & (~(ALIGN_AMT-1)) );
    *((void **)ptr-1) = pa;
    return ptr;
}

inline void aligned_free(void* ptr) {
    if (ptr) std::free(*((void **)ptr-1));
}

template<typename T, int align = ALIGN_AMT>
class Align {
#ifdef USE_ALIGNAS
    alignas(align) T t;
#else
    char buf[sizeof(T) + align - 1];
#endif
public:
    typedef T type;
#ifdef USE_ALIGNAS
    Align() {
    }
    Align(const Align<T>& a) : t(a.t) {
    }
    Align(Align<T>&& a) : t(std::move(a.t)) {
    }
    template<typename... Args>
    Align(Args&&... args) : t(std::forward<Args>(args)...) {
    }

    T* addr() {
        return &t;
    }
    const T* addr() const {
        return &t;
    }
#else
    Align() {
        new (addr()) T();
    }
    Align(const Align<T>& a) {
        new (addr()) T(*a.addr());
    }
    Align(Align<T>&& a) {
        new (addr()) T(std::move(*a.addr()));
    }
    template<typename... Args>
    Align(Args&&... args) {
        new (addr()) T(std::forward<Args>(args)...);
    }
    ~Align() {
        addr()->~T();
    }

    T* addr() {
        return static_cast<T*>(reinterpret_cast<void*>((reinterpret_cast<uintptr_t>(&buf) + align - 1) & ~(align - 1)));
    }
    const T* addr() const {
        return static_cast<const T*>(reinterpret_cast<const void*>((reinterpret_cast<uintptr_t>(&buf) + align - 1) & ~(align - 1)));
    }
#endif

    T* operator->() {
        return addr();
    }
    const T* operator->() const {
        return addr();
    }
    T& operator*() {
        return *addr();
    }
    const T& operator*() const {
        return *addr();
    }

    Align<T>& operator=(const Align<T>& a) {
        *addr() = *a.addr();
        return *this;
    }
    Align<T>& operator=(Align<T>&& a) {
        *addr() = std::move(*a.addr());
        return *this;
    }
};

#endif
