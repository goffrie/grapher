#ifndef _RENDER3D_H_
#define _RENDER3D_H_

#include <initializer_list>
#include <array>
#include <cstring>

#include <QMetaType>
#include <QImage>
#include <qrgb.h>

#include "Expression.h"

#include <boost/config/suffix.hpp>
#include <Vc/Vc>

#include "align.h"

#ifdef ALIGN_TEST
#include <QMessageBox>
#include "stacktrace.h"

template<typename T>
inline void aligncheck(T *_ptr) {
    int dummy = 0;
    quintptr ptr = (quintptr) _ptr;
    if (ptr % alignof(T) != 0) {
        const static QString msg =
                "%1 at address 0x%2, minimum alignment 0x%3\n"
                "Current stack address: 0x%4\n"
                "Expect a crash";
        printStack();
        QMessageBox::critical(0, "Alignment error", msg.arg(
                           QLatin1String(typeid(*_ptr).name()),
                           QString::number(ptr, 16),
                           QString::number(alignof(T), 16),
                           QString::number((quintptr)&dummy, 16)));
    }
}
#endif

template<typename N>
struct Vector3D {
	N v[3];
#ifdef ALIGN_TEST
    Vector3D(N x = 0.f, N y = 0.f, N z = 0.f) {
        aligncheck(this);
        v[0] = x;
        v[1] = y;
        v[2] = z;
    }
#else
    Vector3D(N x = 0.f, N y = 0.f, N z = 0.f):
    	v{x, y, z} { }
#endif
    template<int n> N get() const { return v[n]; }
    template<int n> void set(N w) { v[n] = w; }
    N x() const { return v[0]; }
    N y() const { return v[1]; }
    N z() const { return v[2]; }
    void setX(N x) { v[0] = x; }
    void setY(N y) { v[1] = y; }
    void setZ(N z) { v[2] = z; }
private:
public:
    template<int passes = 2>
    Vector3D<N> normalized() const;
};
Q_DECLARE_METATYPE(Vector3D<Vc::float_v>);
Q_DECLARE_METATYPE(Vector3D<Number>);
static_assert(alignof(Vector3D<Vc::float_v>)%16==0, "Vector3D must be aligned to a 16-byte boundary");

struct Transform3D {
    float rows[4][4];
    Transform3D(const float f[16]) : rows
      {{f[0], f[1], f[2], f[3]},
       {f[4], f[5], f[6], f[7]},
       {f[8], f[9], f[10], f[11]},
       {f[12], f[13], f[14], f[15]}} { }
    Transform3D(std::initializer_list<float> f) : Transform3D(f.begin()) { }
    Transform3D() : rows{{1.f, 0.f, 0.f, 0.f},
                         {0.f, 1.f, 0.f, 0.f},
                         {0.f, 0.f, 1.f, 0.f},
                         {0.f, 0.f, 0.f, 1.f}} { }
    float& operator()(int y, int x) { return rows[y][x]; }
    const float& operator()(int y, int x) const { return rows[y][x]; }
    Transform3D inverted() const;
    Transform3D fit(int w, int h, float x1, float x2, float y1, float y2, float z1, float z2) const;
    static Transform3D translator(float dx, float dy, float dz);
    static Transform3D scaler(float x, float y, float z);
    static Transform3D rotatorX(float theta);
    static Transform3D rotatorY(float theta);
    static Transform3D rotatorZ(float theta);
    const static Transform3D isometricTransform;
};
// static_assert(alignof(Transform3D)%16==0, "Transform3D must be aligned to a 16-byte boundary");
Q_DECLARE_METATYPE(Transform3D);

template <typename N>
inline N dot(const Vector3D<N>& a, const Vector3D<N>& b) {
    return a.x() * b.x()
         + a.y() * b.y()
         + a.z() * b.z();
}

template<typename N>
Vector3D<N> operator*(const Transform3D& t, const Vector3D<N>& v) {
    const N x = v.x(), y = v.y(), z = v.z();
    return Vector3D<N>(
            t(0, 0) * x + t(0, 1) * y + t(0, 2) * z + t(0, 3),
            t(1, 0) * x + t(1, 1) * y + t(1, 2) * z + t(1, 3),
            t(2, 0) * x + t(2, 1) * y + t(2, 2) * z + t(2, 3));
}
template<typename N>
Vector3D<N> operator+(const Vector3D<N>& a, const Vector3D<N>& b) {
    return Vector3D<N>(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
}
template<typename N>
Vector3D<N> operator-(const Vector3D<N>& a, const Vector3D<N>& b) {
    return Vector3D<N>(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
}
template<typename N>
Vector3D<N> operator*(const Vector3D<N>& a, N b) {
    return Vector3D<N>(a.x() * b, a.y() * b, a.z() * b);
}

Transform3D operator*(const Transform3D& a, const Transform3D& b);

bool operator==(const Vector3D<float>& a, const Vector3D<float>& b);
inline bool qFuzzyCompare(const Vector3D<float>& a, const Vector3D<float>& b) { return a==b; }
inline bool operator!=(const Vector3D<float>& a, const Vector3D<float>& b) { return !(a==b); }

QDebug operator<<(QDebug a, const Vector3D<float>& b);
QDebug operator<<(QDebug a, const Vector3D<Vc::float_v>& b);
QDebug operator<<(QDebug a, const Transform3D& b);

class Buffer3D {
    std::size_t m_width, m_height;
    QRgb* m_pixels;
    Vector m_zbuffer;
    float m_colorf[3];
    Transform3D m_transform;
    Vector3D<float> m_viewer, m_light, m_half;
public:
    Buffer3D();
    Buffer3D(std::size_t w, std::size_t h, Transform3D transform);
    Buffer3D(Buffer3D&& buf) { *this = std::move(buf); }
    ~Buffer3D();

    Buffer3D copy() const;

    Buffer3D& operator=(Buffer3D&& buf);

    void clear();

    const Transform3D& transform() const { return m_transform; }

    void drawTransformPoint(const Vector3D<float>& p, QRgb c);
    void drawTransformLine(const Vector3D<float>& p1, const Vector3D<float>& p2, QRgb c);
    void drawTransformPoint(const Vector3D<Vc::float_v>& p, QRgb c);
    void drawTransformLine(const Vector3D<Vc::float_v>& p1, const Vector3D<Vc::float_v>& p2, QRgb c);
    // don't have to be unit vectors
    void setColor(QRgb c);
    void setLight(Vector3D<float> light);
    void drawTransformLitPoint(Vector3D<float> p, Vector3D<float> normal, int idx = -1);
    void drawTransformLitPoint(Vector3D<Vc::float_v> p, Vector3D<Vc::float_v> normal, int idx = -1);


    void setPixel(const Vector3D<float>& p, QRgb c);
    void setPixel(const Vector3D<Vc::float_v>& p, QRgb c);
    void drawLine(const Vector3D<float>& p1, const Vector3D<float>& p2, QRgb c);
    void drawLine(const Vector3D<Vc::float_v>& p1, const Vector3D<Vc::float_v>& p2, QRgb c);
    void drawBuffer(int x, int y, const Buffer3D& buf);

    QImage image() const;
};

// Extracts the first vector out of a vector of vectorized floats.
inline Vector3D<float> vecToScal(const Vector3D<Vc::float_v>& v) {
    return Vector3D<float>(v.v[0][0], v.v[1][0], v.v[2][0]);
}
// Changes floats to float_vs by broadcasting.
inline Vector3D<Vc::float_v> scalToVec(const Vector3D<float>& v) {
    return Vector3D<Vc::float_v>(v.v[0], v.v[1], v.v[2]);
}

template<int passes> void newtonSqrtPass(Vc::float_v& invsqrtnv, Vc::float_v& halfnv) {
    const static Vc::float_v threehalf = 1.5f;
    invsqrtnv *= threehalf - (halfnv * invsqrtnv * invsqrtnv);
    newtonSqrtPass<passes-1>(invsqrtnv, halfnv);
}
template<> inline void newtonSqrtPass<0>(Vc::float_v& invsqrtnv, Vc::float_v& halfnv) {
}

template<int passes>
Vector3D<Vc::float_v> normalize(const Vector3D<Vc::float_v>& v) {
    Vc::float_v sum = v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
    Vc::float_v invsqrt = Vc::rsqrt(sum);
    const static Vc::float_v half = 0.5f;
    Vc::float_v halfsum = half * sum;
    newtonSqrtPass<passes>(invsqrt, halfsum);
    return v * invsqrt;
}
template<int passes> // ignored: always use good precision sqrt anyway
Vector3D<float> normalize(const Vector3D<float>& v) {
    float sum = v.x() * v.x() + v.y() * v.y() + v.z() * v.z();
    float invsqrt = 1.f / sqrt(sum);
    return v * invsqrt;
}
template<typename N> template<int passes> Vector3D<N> Vector3D<N>::normalized() const {
    return normalize<passes>(*this);
}

#endif
