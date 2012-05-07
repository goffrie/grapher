#ifndef _RENDER3D_H_
#define _RENDER3D_H_

#include <initializer_list>
#include <array>

#include <QMatrix4x4>
#include <QMetaType>
#include <QImage>
#include <qrgb.h>

#include "Expression.h"

#include <mmintrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>

#ifdef __SSE3__
#include <pmmintrin.h>
#endif

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

typedef __m128 v4sf;
#define _mm_shufd(xmm, mask) (_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(xmm), mask)))

#ifdef BOOST_NO_CONSTEXPR
#define constexpr
#endif

inline constexpr v4sf vec4(float a, float b, float c, float d) {
    return (v4sf){a, b, c, d};
}
struct Vector3D {
    v4sf v; // { x, y, z, 1 }
    constexpr Vector3D() : v(vec4(0.f, 0.f, 0.f, 1.f)) { }
    constexpr Vector3D(v4sf _v) : v(_v) { }
    constexpr Vector3D(float x, float y, float z) : v(vec4(x, y, z, 1.f)) { }
    template<int n> float get() const { return v[n]; }
    template<int n> void set(float w) { v[n] = w; }
    float x() const { return v[0]; }
    float y() const { return v[1]; }
    float z() const { return v[2]; }
    void setX(float x) { v[0] = x; }
    void setY(float y) { v[1] = y; }
    void setZ(float z) { v[2] = z; }
    //QVector2D toVector2D() const { return QVector2D(v.m[0], v.m[1]); }
    template<int passes = 2>
    Vector3D normalized() const;
    template<int passes>
    static void newtonSqrtPass(v4sf& invsqrtnv, v4sf& halfnv);
};
Q_DECLARE_METATYPE(Vector3D);

struct Transform3D {
    std::array<v4sf, 4> rows;
    constexpr Transform3D(const float f[16]) : rows(
      {{vec4(f[0], f[1], f[2], f[3]),
        vec4(f[4], f[5], f[6], f[7]),
        vec4(f[8], f[9], f[10], f[11]),
        vec4(f[12], f[13], f[14], f[15])}}) { }
    constexpr Transform3D(std::initializer_list<float> f) : Transform3D(f.begin()) { }
    constexpr Transform3D() : rows(
      {{vec4(1.f, 0.f, 0.f, 0.f),
        vec4(0.f, 1.f, 0.f, 0.f),
        vec4(0.f, 0.f, 1.f, 0.f),
        vec4(0.f, 0.f, 0.f, 1.f)}}) { }
    Transform3D inverted() const;
    Transform3D fit(int w, int h, float x1, float x2, float y1, float y2, float z1, float z2) const;
    static Transform3D translator(float dx, float dy, float dz);
    static Transform3D scaler(float x, float y, float z);
    static Transform3D rotatorX(float theta);
    static Transform3D rotatorY(float theta);
    static Transform3D rotatorZ(float theta);
    const static Transform3D isometricTransform;
};
Q_DECLARE_METATYPE(Transform3D);

inline v4sf dot4(const Vector3D& a, const Vector3D& b) {
#ifdef __SSE4_1__
    return _mm_dp_ps(a.v, b.v, 0x7F);
#else
    const static v4sf ones = {1.f, 1.f, 1.f, 1.f};
    __m128 t = _mm_mul_ps(a.v, b.v); // {x1*x2, y1*y2, z1*z2, 1}
    t = _mm_add_ps(t, _mm_shufd(t, 0x4E));
                      // z1*z2, 1, x1*x2, y1*y2
        // x1*x2+z1*z2, y1*y2+1, z1*z2+x1*x2, 1+y1*y2
    return _mm_sub_ps(_mm_add_ps(t, _mm_shufd(t, 0x11)), ones);
                                    // {y1*y2+1, x1*x2+z1*z2, y1*y2+1, x1*x2+z1*z2
                      // {x1*x2+y1*y2+z1*z2+1, x1*x2+y1*y2+z1*z2+1, x1*x2+y1*y2+z1*z2+1, x1*x2+y1*y2+z1*z2+1}
           // {x1*x2+y1*y2+z1*z2} * 4
#endif
}
inline float dot(const Vector3D& a, const Vector3D& b) {
    return _mm_cvtss_f32(dot4(a, b));
}

#ifdef __SSE4_1__
inline Vector3D operator*(const Transform3D& t, const Vector3D& v) {
    return _mm_or_ps(
        _mm_or_ps(
            _mm_dp_ps(t.rows[0], v.v, 0xF1), // store into first field
            _mm_dp_ps(t.rows[1], v.v, 0xF2) // store into second field
        ),
        _mm_or_ps(
            _mm_dp_ps(t.rows[2], v.v, 0xF4), // store into third field
            vec4(0.f, 0.f, 0.f, 1.f)
        )
    );
}
#else
Vector3D operator*(const Transform3D& t, const Vector3D& v);
#endif
Transform3D operator*(const Transform3D& a, const Transform3D& b);
Vector3D operator+(const Vector3D& a, const Vector3D& b);
Vector3D operator-(const Vector3D& a, const Vector3D& b);
Vector3D operator*(const Vector3D& a, float b);

bool operator==(const Vector3D& a, const Vector3D& b);
inline bool qFuzzyCompare(const Vector3D& a, const Vector3D& b) { return a==b; }
inline bool operator!=(const Vector3D& a, const Vector3D& b) { return !(a==b); }

QDebug operator<<(QDebug a, const Vector3D& b);
QDebug operator<<(QDebug a, const Transform3D& b);

class Buffer3D {
    std::size_t m_width, m_height;
    QRgb* m_pixels;
    Vector m_zbuffer;
    Transform3D m_transform;
    Vector3D m_viewer, m_light, m_half;
    v4sf m_colorf;
public:
    Buffer3D();
    Buffer3D(std::size_t w, std::size_t h, Transform3D transform);
    Buffer3D(Buffer3D&& buf) { *this = std::move(buf); }
    ~Buffer3D();

    Buffer3D copy() const;

    Buffer3D& operator=(Buffer3D&& buf);

    void clear();

    const Transform3D& transform() const { return m_transform; }

    void drawTransformPoint(const Vector3D& p, QRgb c);
    void drawTransformLine(const Vector3D& p1, const Vector3D& p2, QRgb c);
    // don't have to be unit vectors
    void setColor(QRgb c);
    void setLight(Vector3D light);
    void drawTransformLitPoint(Vector3D p, Vector3D normal, int idx = -1);


    void setPixel(const Vector3D& p, QRgb c);
    void drawLine(const Vector3D& p1, const Vector3D& p2, QRgb c);
    void drawBuffer(int x, int y, const Buffer3D& buf);

    QImage image() const;
};

template<>
inline void Vector3D::newtonSqrtPass<0>(v4sf& invsqrtnv, v4sf& halfnv) {
}
template<int passes>
void Vector3D::newtonSqrtPass(v4sf& invsqrtnv, v4sf& halfnv) {
    const static v4sf threehalf = {1.5f, 1.5f, 1.5f, 1.5f};
    invsqrtnv *= threehalf - (halfnv * invsqrtnv * invsqrtnv);
    Vector3D::newtonSqrtPass<passes-1>(invsqrtnv, halfnv);
}

#ifdef __SSE3__
template<int passes>
Vector3D Vector3D::normalized() const {
    v4sf vv = v * v; // { x*x, y*y, z*z, 1 }
    vv[3] = 0; // { x*x, y*y, z*z, 0 }
    v4sf hv = _mm_hadd_ps(vv, vv); // {x*x + y*y, z*z + 0} x 2
    v4sf nv = _mm_hadd_ps(hv, hv); // {x*x + y*y + z*z} x 4
    v4sf invsqrtnv = _mm_rsqrt_ps(nv);
    const static v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    v4sf halfnv = half * nv;
    newtonSqrtPass<passes>(invsqrtnv, halfnv);
    Vector3D ret(v * invsqrtnv);
    ret.v[3] = 1.f;
    return ret;
}
#else
template<int passes>
Vector3D Vector3D::normalized() const {
    v4sf nv = dot4(v, v); // {x*x + y*y + z*z} x 4
    v4sf invsqrtnv = _mm_rsqrt_ps(nv);
    const static v4sf half = {0.5f, 0.5f, 0.5f, 0.5f};
    v4sf halfnv = half * nv;
    newtonSqrtPass<passes>(invsqrtnv, halfnv);
    Vector3D ret(v * invsqrtnv);
    ret.v[3] = 1.f;
    return ret;
}
#endif


#endif
