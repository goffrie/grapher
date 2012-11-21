#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <boost/config.hpp>
#include <cstdint>
#include <Vc/Vc>
#include <align.h>

typedef float Number;
typedef std::size_t uz;
typedef std::uint64_t u64;
typedef std::int64_t i64;
typedef std::uint32_t u32;
typedef std::int32_t i32;
typedef std::uint8_t u8;

BOOST_CONSTEXPR_OR_CONST uint SSE_VECTOR_SIZE = Vc::float_v::Size;

#define VECTOR_ALLOC(num) ((Vector) aligned_malloc((((num) + SSE_VECTOR_SIZE - 1) & (~(Vc::float_v::Size - 1)))*sizeof(Number)))
#define VECTOR_FREE(ptr) aligned_free(ptr)

#define VECTOR_LOOP for (uint i = 0; i < size; i += Vc::float_v::Size)
#define V(a) (*reinterpret_cast<Vc::float_v*>(a+i))

typedef Number __attribute__((aligned(32))) * Vector;
typedef Number __attribute__((aligned(32))) * __restrict VectorR;

struct VectorDeleter {
    void operator()(Vector ptr) const throw() {
        VECTOR_FREE(ptr);
    }
};

typedef std::unique_ptr<Number[], VectorDeleter> UVector;

static_assert(sizeof(UVector) == sizeof(Vector), "UVector is not the right size!");

#endif
