#include "Render3D.h"

#include <QVector2D>
#include <QDebug>

#include <cstring>

#include <gsl/gsl_math.h>

#include <boost/config/suffix.hpp>

BOOST_CONSTEXPR_OR_CONST
float rt2 = M_SQRT2,
      rt3 = M_SQRT3;

const Transform3D Transform3D::isometricTransform = Transform3D{
    1.f/rt2, 0.f, -1.f/rt2, 0.f,
    1.f/(rt2*rt3), rt2/rt3, 1.f/(rt2*rt3), 0.f,
    1.f/rt3, -1.f/rt3, 1.f/rt3, 0.f,
    0.f, 0.f, 0.f, 1.f
};

Transform3D Transform3D::fit(int w, int h, float x1, float x2, float y1, float y2, float z1, float z2) const {
    const Vector3D<float> points[8] = {
        Vector3D<float>(x1, y1, z1),
        Vector3D<float>(x1, y1, z2),
        Vector3D<float>(x1, y2, z1),
        Vector3D<float>(x1, y2, z2),
        Vector3D<float>(x2, y1, z1),
        Vector3D<float>(x2, y1, z2),
        Vector3D<float>(x2, y2, z1),
        Vector3D<float>(x2, y2, z2)
    };
    float left, right, top, bottom;
    for (int i = 0; i < 8; ++i) {
        Vector3D<float> p = (*this) * points[i];
        if (i == 0 || left > p.x()) left = p.x();
        if (i == 0 || right < p.x()) right = p.x();
        if (i == 0 || top > p.y()) top = p.y();
        if (i == 0 || bottom < p.y()) bottom = p.y();
    }
    float scale = qMin(w / (right - left), h / (bottom - top));
    float dx = (w - ((right - left) * scale)) / 2 - left * scale,
          dy = (h - ((bottom - top) * scale)) / 2 - top * scale;
    return Transform3D::translator(dx, dy, 0.) * Transform3D::scaler(scale, scale, 1) * (*this);
}

Transform3D Transform3D::translator(float dx, float dy, float dz) {
    Transform3D m;
    m.rows[0][3] = dx;
    m.rows[1][3] = dy;
    m.rows[2][3] = dz;
    return m;
}

Transform3D Transform3D::scaler(float dx, float dy, float dz) {
    Transform3D m;
    m.rows[0][0] = dx;
    m.rows[1][1] = dx;
    m.rows[2][2] = dx;
    return m;
}

Transform3D Transform3D::rotatorX(float theta) {
    const float c = std::cos(theta), s = std::sin(theta);
    return Transform3D{
        1, 0, 0, 0,
        0, c,-s, 0,
        0, s, c, 0,
        0, 0, 0, 1};
}

Transform3D Transform3D::rotatorY(float theta) {
    const float c = std::cos(theta), s = std::sin(theta);
    return Transform3D{
        c, 0, s, 0,
        0, 1, 0, 0,
       -s, 0, c, 0,
        0, 0, 0, 1};
}

Transform3D Transform3D::rotatorZ(float theta) {
    const float c = std::cos(theta), s = std::sin(theta);
    return Transform3D{
        c,-s, 0, 0,
        s, c, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1};
}

// matrix multiplication
Transform3D operator*(const Transform3D& a, const Transform3D& b) {
    Transform3D result{0.f, 0.f, 0.f, 0.f,
                       0.f, 0.f, 0.f, 0.f,
                       0.f, 0.f, 0.f, 0.f,
                       0.f, 0.f, 0.f, 0.f};
    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            for (int i = 0; i < 4; ++i) {
                result(y, x) += a(y, i) * b(i, x);
            }
        }
    }
    return result;
}

bool operator==(const Vector3D<float>& a, const Vector3D<float>& b) {
    return qFuzzyCompare(a.v[0], b.v[0]) && qFuzzyCompare(a.v[1], b.v[1]) && qFuzzyCompare(a.v[2], b.v[2]);
}

QDebug operator<<(QDebug a, const Vector3D<float>& b) {
    a.nospace() << "Vector3D(" << b.x() << ',' << b.y() << ',' << b.z() << ')';
    return a.space();
}
QDebug operator<<(QDebug a, const Vector3D<Vc::float_v>& b) {
    auto q = a.nospace() << "Vector3D<float_v>{";
    for (int i = 0; i < 4; ++i) {
        q = q << '(' << b.x()[i] << ',' << b.y()[i] << ',' << b.z()[i] << ((i==3)?"),":")");
    }
    return a.space();
}
QDebug operator<<(QDebug a, const Transform3D& b) {
    a.nospace() << "Transform3D{\n"
        << b.rows[0][0] << ',' << b.rows[0][1] << ',' << b.rows[0][2] << ',' << b.rows[0][3] << '\n'
        << b.rows[1][0] << ',' << b.rows[1][1] << ',' << b.rows[1][2] << ',' << b.rows[1][3] << '\n'
        << b.rows[2][0] << ',' << b.rows[2][1] << ',' << b.rows[2][2] << ',' << b.rows[2][3] << '\n'
        << b.rows[3][0] << ',' << b.rows[3][1] << ',' << b.rows[3][2] << ',' << b.rows[3][3] << '\n'
        << '}';
    return a.space();
}


Buffer3D::Buffer3D() : m_width(0), m_height(0), m_pixels(NULL), m_zbuffer(NULL) {
}

Buffer3D::Buffer3D(std::size_t w, std::size_t h, Transform3D transform) : m_width(w), m_height(h),
m_pixels((QRgb*) aligned_malloc(((w*h + (Vc::float_v::Size-1)) & ~(Vc::float_v::Size-1))*sizeof(QRgb))), m_zbuffer(VECTOR_ALLOC(w*h)), m_transform(transform) {
    clear();
    Transform3D inv = transform.inverted();
    m_viewer = Vector3D<float>(inv(0, 2), inv(1, 2), inv(2, 2)).normalized<6>();
}

Buffer3D::~Buffer3D() {
    aligned_free(m_pixels);
    VECTOR_FREE(m_zbuffer);
}

Buffer3D& Buffer3D::operator=(Buffer3D&& buf) {
    m_width = buf.m_width;
    m_height = buf.m_height;
    std::swap(m_pixels, buf.m_pixels);
    std::swap(m_zbuffer, buf.m_zbuffer);
    memcpy(m_colorf, buf.m_colorf, sizeof(m_colorf));
    m_transform = buf.m_transform;
    m_viewer = buf.m_viewer;
    m_light = buf.m_light;
    m_half = buf.m_half;
    return *this;
}

Buffer3D Buffer3D::copy() const {
    Buffer3D r;
    r.m_width = m_width;
    r.m_height = m_height;
    memcpy(r.m_colorf, m_colorf, sizeof(m_colorf));
    r.m_transform = m_transform;
    r.m_viewer = m_viewer;
    r.m_light = m_light;
    r.m_half = m_half;
    const int size = m_width*m_height;
    r.m_pixels = (QRgb*) aligned_malloc(((size + (Vc::float_v::Size-1)) & ~(Vc::float_v::Size-1))*sizeof(QRgb));
    r.m_zbuffer = VECTOR_ALLOC(size);
    std::memcpy(r.m_pixels, m_pixels, sizeof(QRgb)*size);
    std::memcpy(r.m_zbuffer, m_zbuffer, sizeof(Number)*size);
    return r;
}

void Buffer3D::clear() {
    std::memset(m_pixels, -1, m_width*m_height*sizeof(QRgb));
    uz size = m_width * m_height;
    const Vc::float_v neginf(GSL_NEGINF);
    for (uz i = 0; i < size; i += Vc::float_v::Size) {
        neginf.store(m_zbuffer+i);
    }
}

void Buffer3D::drawTransformPoint(const Vector3D<float>& p, QRgb c) {
    const Vector3D<float> tp = m_transform * p;
    const int x = qRound(tp.x()), y = qRound(tp.y());
    const Number z = tp.z();
    if (x < 0 || x >= m_width || y < 0 || y >= m_height) return;
    const std::size_t idx = y * m_width + x;
    if (z < m_zbuffer[idx]) return;
    m_zbuffer[idx] = z;
    m_pixels[idx] = c;
}

void Buffer3D::setPixel(const Vector3D<float>& p, QRgb c) {
    const int x = qRound(p.x()), y = qRound(p.y());
    const Number z = p.z();
    if (x < 0 || x >= m_width || y < 0 || y >= m_height) return;
    const std::size_t idx = y * m_width + x;
    if (z < m_zbuffer[idx]) return;
    m_zbuffer[idx] = z;
    m_pixels[idx] = c | 0xff000000;
}

void Buffer3D::drawTransformLine(const Vector3D<float>& p1, const Vector3D<float>& p2, QRgb c) {
    drawLine(m_transform * p1, m_transform * p2, c);
}

void Buffer3D::drawLine(const Vector3D<float>& p1, const Vector3D<float>& p2, QRgb c) {
    Vector3D<float> dp = p2 - p1;
    float tstep = 1.f / qMax(qAbs(dp.x()), qAbs(dp.y()));
    for (float t = 0; t <= 1; t += tstep) {
        setPixel(p1 + dp*t, c);
    }
}

void Buffer3D::setColor(QRgb c) {
    m_colorf[0] = qRed(c);
    m_colorf[1] = qGreen(c);
    m_colorf[2] = qBlue(c);
}

void Buffer3D::setLight(Vector3D<float> light) {
    m_light = light.normalized<6>();
    m_half = (m_light + m_viewer).normalized<6>();
}
/*
inline
QDebug
operator<<(QDebug d, const v4sf& f) {
    d.nospace() << '(' << f[0] << ',' << f[1] << ',' << f[2] << ',' << f[3] << ')';
    return d.space();
}
*/
float rsqrt(float n) {
    return 1.f / sqrt(n);
}
unsigned char clamp(float n) {
    return (n < 0) ? 0 : (n > 255) ? 255 : n;
}

void Buffer3D::drawTransformLitPoint(Vector3D<float> p, Vector3D<float> n, int idx) {
    const Vector3D<float> tp = m_transform * p;
    const float z = tp.z();
    if (idx == -1) {
        const int x = qRound(tp.x()), y = qRound(tp.y());
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) return;
        idx = y * m_width + x;
    }
    if (z < m_zbuffer[idx]) return;
    m_zbuffer[idx] = z;
    float rsqrt_nn = rsqrt(dot(n, n)),
          nv = dot(n, m_viewer),
          nl = dot(n, m_light),
          nh = dot(n, m_half);
    if (nv < 0) {
        // If `n' is facing the wrong way relative to the viewer, then flip everything.
        nl = -nl;
        nh = -nh;
    }

    float lighting = 0.2f; // ambient lighting
    lighting += std::max(rsqrt_nn * nl, 0.f); // diffuse term = n * l / |n||l| = n * l / |n| = cos(theta)
    float litc[3] = { lighting * m_colorf[0], lighting * m_colorf[1], lighting * m_colorf[2] };
    float rvd = std::max(rsqrt_nn * nh, 0.f); // n * h / |n|
    BOOST_STATIC_CONSTEXPR float specbase[3] = {150.f, 150.f, 150.f};
    float rvd2 = rvd * rvd, rvd4 = rvd2 * rvd2, rvd8 = rvd4 * rvd4;
    litc[0] += specbase[0] * rvd8;
    litc[1] += specbase[1] * rvd8;
    litc[2] += specbase[2] * rvd8;
    m_pixels[idx] = qRgb(clamp(litc[0]), clamp(litc[1]), clamp(litc[2]));
}

void Buffer3D::drawBuffer(int x, int y, const Buffer3D& buf) {
    int ymax = qMin(y + buf.m_height, m_height);
    int xmax = qMin(x + buf.m_width, m_width);
    for (int Y = y; Y < ymax; ++Y) {
        for (int X = x; X < xmax; ++X) {
            int idx = Y*m_width + X;
            bool copy = m_zbuffer[idx] < buf.m_zbuffer[idx];
            m_zbuffer[idx] = copy ? buf.m_zbuffer[idx] : m_zbuffer[idx];
            m_pixels[idx] = copy ? buf.m_pixels[idx] : m_pixels[idx];
        }
    }
}

QImage Buffer3D::image() const {
    return QImage(reinterpret_cast<const uchar*>(m_pixels), (int)m_width, (int)m_height, (int)(m_width*sizeof(QRgb)), QImage::Format_RGB32);
}
