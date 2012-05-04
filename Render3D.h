#ifndef _RENDER3D_H_
#define _RENDER3D_H_

#include <initializer_list>

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

typedef __m128 v4sf;
#define _mm_shufd(xmm, mask) (_mm_castsi128_ps(_mm_shuffle_epi32(_mm_castps_si128(xmm), mask)))

inline v4sf vec4(float a, float b, float c, float d) {
    const v4sf r = {a, b, c, d};
    return r;
}
struct alignas(16) Vector3D {
    v4sf v; // { x, y, z, 1 }
    Vector3D() : v(vec4(0.f, 0.f, 0.f, 1.f)) { }
    Vector3D(v4sf _v) : v(_v) { }
    Vector3D(float x, float y, float z) : v(vec4(x, y, z, 1.f)) { }
    template<int n> float get() const { return v[n]; }
    template<int n> void set(float w) const { v[n] = w; }
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
struct alignas(16) Transform3D {
    v4sf rows[4];
    Transform3D(const float* f);
    Transform3D(const std::initializer_list<float>& f = {1.f, 0.f, 0.f, 0.f,
        0.f, 1.f, 0.f, 0.f,
        0.f, 0.f, 1.f, 0.f,
        0.f, 0.f, 0.f, 1.f}) : Transform3D(f.begin()) { }
    Transform3D inverted(bool* invertible) const;
    static Transform3D translator(qreal dx, qreal dy, qreal dz);
    static Transform3D scaler(qreal x, qreal y, qreal z);
};
Q_DECLARE_METATYPE(Transform3D);

inline v4sf dot4(const Vector3D& a, const Vector3D& b) {
    const static v4sf ones = {1.f, 1.f, 1.f, 1.f};
    __m128 t = _mm_mul_ps(a.v, b.v); // {x1*x2, y1*y2, z1*z2, 1}
    t = _mm_add_ps(t, _mm_shufd(t, 0x4E));
                      // z1*z2, 1, x1*x2, y1*y2
        // x1*x2+z1*z2, y1*y2+1, z1*z2+x1*x2, 1+y1*y2
    return _mm_sub_ps(_mm_add_ps(t, _mm_shufd(t, 0x11)), ones);
                                    // {y1*y2+1, x1*x2+z1*z2, y1*y2+1, x1*x2+z1*z2
                      // {x1*x2+y1*y2+z1*z2+1, x1*x2+y1*y2+z1*z2+1, x1*x2+y1*y2+z1*z2+1, x1*x2+y1*y2+z1*z2+1}
           // {x1*x2+y1*y2+z1*z2} * 4
}
inline float dot(const Vector3D& a, const Vector3D& b) {
    return _mm_cvtss_f32(dot4(a, b));
}

Vector3D operator*(const Transform3D& t, const Vector3D& v);
Transform3D operator*(const Transform3D& a, const Transform3D& b);
Vector3D operator+(const Vector3D& a, const Vector3D& b);
Vector3D operator-(const Vector3D& a, const Vector3D& b);
Vector3D operator*(const Vector3D& a, float b);

bool operator==(const Vector3D& a, const Vector3D& b);
inline bool qFuzzyCompare(const Vector3D& a, const Vector3D& b) { return a==b; }
inline bool operator!=(const Vector3D& a, const Vector3D& b) { return !(a==b); }

QDebug operator<<(QDebug a, const Vector3D& b);
QDebug operator<<(QDebug a, const Transform3D& b);

Transform3D isometricTransform(int w, int h, int x1, int x2, int y1, int y2, int z1, int z2);

class Buffer3D {
    std::size_t m_width, m_height;
    QRgb* m_pixels;
    Number* m_zbuffer;
    Transform3D m_transform;
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
    void drawTransformLitPoint(const Vector3D& p, QRgb c, const Vector3D& normal, const Vector3D& light, const Vector3D& eyeRay, int idx = -1);

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
    v4sf vv = { v * v }; // { x*x, y*y, z*z, 1 }
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
