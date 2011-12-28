#ifndef _RENDER3D_H_
#define _RENDER3D_H_

#include <initializer_list>

#include <QMatrix4x4>
#include <QMetaType>
#include <QImage>
#include <qrgb.h>

#include "Expression.h"

#include <xmmintrin.h>

typedef __m128 v4sf;
union v4sfi {
    v4sf v;
    float m[4];
};
inline v4sf vec4(float a, float b, float c, float d) {
    const v4sf r = {a, b, c, d};
    return r;
}
struct Vector3D {
    v4sfi v; // { x, y, z, 1 }
    Vector3D() : v({vec4(0.f, 0.f, 0.f, 1.f)}) { }
    Vector3D(v4sf _v) : v({_v}) { }
    Vector3D(float x, float y, float z) : v({vec4(x, y, z, 1.f)}) { }
    float x() const { return v.m[0]; }
    float y() const { return v.m[1]; }
    float z() const { return v.m[2]; }
    void setX(float x) { v.m[0] = x; }
    void setY(float y) { v.m[1] = y; }
    void setZ(float z) { v.m[2] = z; }
    //QVector2D toVector2D() const { return QVector2D(v.m[0], v.m[1]); }
    Vector3D normalized() const;
};
Q_DECLARE_METATYPE(Vector3D);
struct Transform3D {
    v4sfi rows[4];
    Transform3D(const float* f = identity);
    const static float identity[16];
    Transform3D inverted(bool* invertible) const;
    static Transform3D translator(qreal dx, qreal dy, qreal dz);
    static Transform3D scaler(qreal x, qreal y, qreal z);
};
Q_DECLARE_METATYPE(Transform3D);

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

    Buffer3D& operator=(Buffer3D&& buf);

    void clear();

    const Transform3D& transform() const { return m_transform; }

    void drawPoint(const Vector3D& p, QRgb c);
    // don't have to be unit vectors
    void drawLitPoint(const Vector3D& p, QRgb c, const Vector3D& normal, const Vector3D& l, int idx = -1);

    void setPixel(const Vector3D& p, QRgb c);

    QImage image() const { return QImage(reinterpret_cast<const uchar*>(m_pixels), (int)m_width, (int)m_height, (int)(m_width*sizeof(QRgb)), QImage::Format_RGB32); }
};


#endif
