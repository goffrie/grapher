#include "Render3D.h"

#include <QVector2D>

#include <cstring>

#include <gsl/gsl_math.h>

Transform3D Transform3D::inverted() const {
    __m128 minor0, minor1, minor2, minor3;
    __m128 row0, row1, row2, row3;
    __m128 det, tmp1;
    const float* const src = (const float*) rows.begin();

    tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src)), (__m64*)(src+ 4));
    row1 = _mm_loadh_pi(_mm_loadl_pi(row1, (__m64*)(src+8)), (__m64*)(src+12));
    row0 = _mm_shuffle_ps(tmp1, row1, 0x88);
    row1 = _mm_shuffle_ps(row1, tmp1, 0xDD);
    tmp1 = _mm_loadh_pi(_mm_loadl_pi(tmp1, (__m64*)(src+ 2)), (__m64*)(src+ 6));
    row3 = _mm_loadh_pi(_mm_loadl_pi(row3, (__m64*)(src+10)), (__m64*)(src+14));
    row2 = _mm_shuffle_ps(tmp1, row3, 0x88);
    row3 = _mm_shuffle_ps(row3, tmp1, 0xDD);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row2, row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    minor0 = _mm_mul_ps(row1, tmp1);
    minor1 = _mm_mul_ps(row0, tmp1);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    minor0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), minor0);
    minor1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor1);
    minor1 = _mm_shuffle_ps(minor1, minor1, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row1, row2);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    minor0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor0);
    minor3 = _mm_mul_ps(row0, tmp1);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row3, tmp1));
    minor3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor3);
    minor3 = _mm_shuffle_ps(minor3, minor3, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
    row2 = _mm_shuffle_ps(row2, row2, 0x4E);

    minor0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor0);
    minor2 = _mm_mul_ps(row0, tmp1);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    minor0 = _mm_sub_ps(minor0, _mm_mul_ps(row2, tmp1));
    minor2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), minor2);
    minor2 = _mm_shuffle_ps(minor2, minor2, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row1);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    minor2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor2);
    minor3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), minor3);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    minor2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), minor2);
    minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row2, tmp1));
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row2, tmp1));
    minor2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor2);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    minor1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), minor1);
    minor2 = _mm_sub_ps(minor2, _mm_mul_ps(row1, tmp1));
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row2);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    minor1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), minor1);
    minor3 = _mm_sub_ps(minor3, _mm_mul_ps(row1, tmp1));

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    minor1 = _mm_sub_ps(minor1, _mm_mul_ps(row3, tmp1));
    minor3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), minor3);
    // -----------------------------------------------
    det = _mm_mul_ps(row0, minor0);
    det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
    det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
    tmp1 = _mm_rcp_ss(det);

    det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
    det = _mm_shuffle_ps(det, det, 0x00);

    Transform3D ret;
    float* dest = (float*)ret.rows.begin();

    minor0  = _mm_mul_ps(det, minor0);
    _mm_storel_pi((__m64*)(dest), minor0);
    _mm_storeh_pi((__m64*)(dest+2), minor0);

    minor1  = _mm_mul_ps(det, minor1);
    _mm_storel_pi((__m64*)(dest+4), minor1);
    _mm_storeh_pi((__m64*)(dest+6), minor1);

    minor2  = _mm_mul_ps(det, minor2);
    _mm_storel_pi((__m64*)(dest+ 8), minor2);
    _mm_storeh_pi((__m64*)(dest+10), minor2);

    minor3  = _mm_mul_ps(det, minor3);
    _mm_storel_pi((__m64*)(dest+12), minor3);
    _mm_storeh_pi((__m64*)(dest+14), minor3);

    return ret;
}

const Transform3D Transform3D::isometricTransform{
    1./M_SQRT2, 0., -1./M_SQRT2, 0.,
    1./(M_SQRT2*M_SQRT3), M_SQRT2/M_SQRT3, 1./(M_SQRT2*M_SQRT3), 0.,
    1./M_SQRT3, -1./M_SQRT3, 1./M_SQRT3, 0.,
    0., 0., 0., 1.
};

Transform3D Transform3D::fit(int w, int h, float x1, float x2, float y1, float y2, float z1, float z2) const {
    const Vector3D points[8] = {
        Vector3D(x1, y1, z1),
        Vector3D(x1, y1, z2),
        Vector3D(x1, y2, z1),
        Vector3D(x1, y2, z2),
        Vector3D(x2, y1, z1),
        Vector3D(x2, y1, z2),
        Vector3D(x2, y2, z1),
        Vector3D(x2, y2, z2)
    };
    float left, right, top, bottom;
    for (int i = 0; i < 8; ++i) {
        Vector3D p = (*this) * points[i];
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
    return Transform3D({
        1, 0, 0, 0,
        0, c,-s, 0,
        0, s, c, 0,
        0, 0, 0, 1});
}

Transform3D Transform3D::rotatorY(float theta) {
    const float c = std::cos(theta), s = std::sin(theta);
    return Transform3D({
        c, 0, s, 0,
        0, 1, 0, 0,
       -s, 0, c, 0,
        0, 0, 0, 1});
}

Transform3D Transform3D::rotatorZ(float theta) {
    const float c = std::cos(theta), s = std::sin(theta);
    return Transform3D({
        c,-s, 0, 0,
        s, c, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1});
}

#ifdef __SSE3__
inline __m128 __attribute__((__gnu_inline__, __always_inline__))
haddps(__m128 a, __m128 b) { return _mm_hadd_ps(a, b); }
#else
inline __m128 __attribute__((__gnu_inline__, __always_inline__))
haddps(__m128 a, __m128 b) { return _mm_shuffle_ps(a, b, 0x88) + _mm_shuffle_ps(a, b, 0xDD); }
#endif
#ifndef __SSE4_1__
Vector3D operator*(const Transform3D& t, const Vector3D& v) {
    const v4sf a = t.rows[0] * v.v; // {t00 * v0, t01 * v1, t02 * v2, t03 * v3}
    const v4sf b = t.rows[1] * v.v; // {t10 * v0, t11 * v1, t12 * v2, t13 * v3}
    const v4sf c = t.rows[2] * v.v; // {t20 * v0, t21 * v1, t22 * v2, t23 * v3}
    const v4sf ab = haddps(a, b); // {t00 * v0 + t01 * v1, t02 * v2 + t03 * v3, t10 * v0 + t11 * v1, t12 * v2 + t13 * v3}
    const static v4sf onezeros = {1.f, 0.f, 0.f, 0.f};
    const v4sf c1 = haddps(c, onezeros); // {t20 * v0 + t21 * v1, t22 * v2 + t23 * v3, 1.f, 0.f}
    const v4sf result = haddps(ab, c1); // {t00 * v0 + t01 * v1 + t02 * v2 + t03 * v3,
                                             //  t10 * v0 + t11 * v1 + t12 * v2 + t13 * v3,
                                             //  t20 * v0 + t21 * v1 + t22 * v2 + t23 * v3,
                                             //  1.f}
    return Vector3D(result);
}
#endif

// TODO: optimize????
Transform3D operator*(const Transform3D& a, const Transform3D& b) {
    constexpr static float zeroes[16] = { 0.f };
    Transform3D result(zeroes);
    for (int y = 0; y < 4; ++y) {
        for (int x = 0; x < 4; ++x) {
            for (int i = 0; i < 4; ++i) {
                result.rows[y][x] += a.rows[y][i] * b.rows[i][x];
            }
        }
    }
    return result;
}

Vector3D operator+(const Vector3D& a, const Vector3D& b) {
    Vector3D result(a.v + b.v);
    result.v[3] = 1.f;
    return result;
}
Vector3D operator-(const Vector3D& a, const Vector3D& b) {
    Vector3D result(a.v - b.v);
    result.v[3] = 1.f;
    return result;
}
Vector3D operator*(const Vector3D& a, float b) {
    v4sf B = {b, b, b, 1.f};
    return Vector3D(a.v * B);
}

bool operator==(const Vector3D& a, const Vector3D& b) {
    return qFuzzyCompare(a.v[0], b.v[0]) && qFuzzyCompare(a.v[1], b.v[1]) && qFuzzyCompare(a.v[2], b.v[2]);
}

QDebug operator<<(QDebug a, const Vector3D& b) {
    a.nospace() << "Vector3D(" << b.x() << ',' << b.y() << ',' << b.z() << ")[" << b.v[3] << ']';
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
m_pixels((QRgb*) aligned_malloc(((w*h + 3) & ~3)*sizeof(QRgb))), m_zbuffer(VECTOR_ALLOC(w*h)), m_transform(transform) {
    clear();
    const v4sf dz = {0.f, 0.f, 1.f, 0.f};
    m_viewer = (m_transform.inverted() * dz).normalized<6>();
}

Buffer3D::~Buffer3D() {
    aligned_free(m_pixels);
    VECTOR_FREE(m_zbuffer);
}

Buffer3D& Buffer3D::operator=(Buffer3D&& buf) {
    std::swap(m_width, buf.m_width);
    std::swap(m_height, buf.m_height);
    std::swap(m_pixels, buf.m_pixels);
    std::swap(m_zbuffer, buf.m_zbuffer);
    std::swap(m_transform, buf.m_transform);
    std::swap(m_viewer, buf.m_viewer);
    std::swap(m_light, buf.m_light);
    std::swap(m_half, buf.m_half);
    std::swap(m_colorf, buf.m_colorf);
    return *this;
}

Buffer3D Buffer3D::copy() const {
    Buffer3D r;
    r.m_width = m_width;
    r.m_height = m_height;
    r.m_transform = m_transform;
    r.m_viewer = m_viewer;
    r.m_light = m_light;
    r.m_half = m_half;
    const int size = m_width*m_height;
    r.m_pixels = (QRgb*) aligned_malloc(((size + 3) & ~3)*sizeof(QRgb));
    r.m_zbuffer = VECTOR_ALLOC(size);
    std::memcpy(r.m_pixels, m_pixels, sizeof(QRgb)*size);
    std::memcpy(r.m_zbuffer, m_zbuffer, sizeof(Number)*size);
    return r;
}

void Buffer3D::clear() {
    std::size_t size = (m_width * m_height + 3) / 4;
    constexpr __v4si trans = {-1, -1, -1, -1};
    for (__v4si* pixels = reinterpret_cast<__v4si*>(m_pixels), * end = pixels + size;
        pixels != end;
        ++pixels) *pixels = trans;
    const v4sf neginf = _mm_set_ps1(GSL_NEGINF);
    for (v4sf* zbuffer = reinterpret_cast<v4sf*>(m_zbuffer), * end = zbuffer + size;
        zbuffer != end;
        ++zbuffer) *zbuffer = neginf;
}

void Buffer3D::drawTransformPoint(const Vector3D& p, QRgb c) {
    const Vector3D tp = m_transform * p;
    const int x = qRound(tp.x()), y = qRound(tp.y());
    const Number z = tp.z();
    if (x < 0 || x >= m_width || y < 0 || y >= m_height) return;
    const std::size_t idx = y * m_width + x;
    if (z < m_zbuffer[idx]) return;
    m_zbuffer[idx] = z;
    m_pixels[idx] = c;
}

void Buffer3D::setPixel(const Vector3D& p, QRgb c) {
    const int x = qRound(p.x()), y = qRound(p.y());
    const Number z = p.z();
    if (x < 0 || x >= m_width || y < 0 || y >= m_height) return;
    const std::size_t idx = y * m_width + x;
    if (z < m_zbuffer[idx]) return;
    m_zbuffer[idx] = z;
    m_pixels[idx] = c | 0xff000000;
}

void Buffer3D::drawTransformLine(const Vector3D& p1, const Vector3D& p2, QRgb c) {
    drawLine(m_transform * p1, m_transform * p2, c);
}

void Buffer3D::drawLine(const Vector3D& p1, const Vector3D& p2, QRgb c) {
    Vector3D dp = p2 - p1;
    float tstep = 1.f / qMax(qAbs(dp.x()), qAbs(dp.y()));
    for (float t = 0; t <= 1; t += tstep) {
        setPixel(p1 + dp*t, c);
    }
}

void Buffer3D::setColor(QRgb c) {
    m_colorf = _mm_cvtpu8_ps(_mm_cvtsi32_si64(c));
}

void Buffer3D::setLight(Vector3D light) {
    m_light = light.normalized<6>();
    m_half = (m_light + m_viewer).normalized<6>();
    m_light.v[3] = 0.f;
    m_half.v[3] = 0.f;
}

QDebug
operator<<(QDebug d, v4sf f) {
    d.nospace() << '(' << f[0] << ',' << f[1] << ',' << f[2] << ',' << f[3] << ')';
    return d.space();
}

void Buffer3D::drawTransformLitPoint(Vector3D p, Vector3D normal, int idx) {
    const Vector3D tp = m_transform * p;
    const Number z = tp.z();
    if (idx == -1) {
        const int x = qRound(tp.x()), y = qRound(tp.y());
        if (x < 0 || x >= m_width || y < 0 || y >= m_height) return;
        idx = y * m_width + x;
    }
    if (z < m_zbuffer[idx]) return;
    m_zbuffer[idx] = z;
    static const __v4si nmask = {-1,-1,-1,0};
    static const __v4si invneg = {-1<<31,-1<<31,-1<<31,-1<<31};
    static const v4sf zero = {0.f,0.f,0.f,0.f};
    v4sf n = _mm_and_si128(normal.v, nmask);
    const v4sf _nlvhd = haddps(haddps(n * n, n * m_light.v), haddps(n * m_viewer.v, n * m_half.v)); // { n * n, n * l, n * v, n * h }
    const v4sf rcps_nlvhd = _mm_rsqrt_ps(_nlvhd); // { 1/|n|, 1/sqrt(n*l), 1/sqrt(n*v), 1/sqrt(n*h) }
    const v4sf nlvhd = _mm_xor_ps(_nlvhd, _mm_and_ps(_mm_cmplt_ps(_mm_shuffle_ps(_nlvhd, _nlvhd, 0xAA), zero), invneg)); // invert if n*v < 0
    float lighting = 0.2f; // ambient lighting
    lighting += std::max(rcps_nlvhd[0] * nlvhd[1], 0.f); // diffuse term = n * l / |n||l| = n * l / |n| = cos(theta)
    const v4sf vlighting = _mm_set_ps1(lighting);
    v4sf litc = vlighting * m_colorf;
    const v4sf rvd = _mm_max_ps(_mm_set_ps1(rcps_nlvhd[0] * nlvhd[3]), zero); // n * h / |n|
    static const v4sf specbase = {150.f, 150.f, 150.f, 0.f};
    const v4sf rvd2 = rvd * rvd;
    const v4sf rvd4 = rvd2 * rvd2;
    litc += specbase * rvd4 * rvd4;
    const __v4si litcw = _mm_cvttps_epi32(litc);
    const __m128i litcw16 = _mm_packs_epi32(litcw, litcw);
    const __v4si litcw8 = _mm_packus_epi16(litcw16, litcw16);
    m_pixels[idx] = litcw8[0] | 0xff000000;
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
