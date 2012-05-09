#ifndef ALIGN_H
#define ALIGN_H

#include <cstdlib>
#include <cstdint>

inline void* aligned_malloc(size_t size) {
    void* pa = std::malloc(size + 15 + sizeof(void *));
    if (!pa) return NULL;
    void* ptr = (void*) ( ((uintptr_t)pa + sizeof(void *) + 15) & (~15) );
    *((void **)ptr-1) = pa;
    return ptr;
}

inline void aligned_free(void* ptr) {
    if (ptr) std::free(*((void **)ptr-1));
}

template<typename T, int align = 16>
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
        return reinterpret_cast<T*>((reinterpret_cast<uintptr_t>(&buf) + align - 1) & ~(align - 1));
    }
    const T* addr() const {
        return reinterpret_cast<const T*>((reinterpret_cast<uintptr_t>(&buf) + align - 1) & ~(align - 1));
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
