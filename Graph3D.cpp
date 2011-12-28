#include "Graph3D.h"

#include <QtConcurrentRun>
#include <QtConcurrentMap>
#include <QThreadStorage>

#include <iostream>
#include <functional>

#include <gsl/gsl_sys.h>
#include <gsl/gsl_nan.h>

Graph3D::Graph3D(QObject* parent) : Graph(parent) {

}

void Graph3D::setupRestart(const Transform3D& t, int width, int height, Vector3D boxa, Vector3D boxb, Vector3D light) {
    cancel();
    m_width = width;
    m_height = height;
    m_boxa = boxa;
    m_boxb = boxb;
    m_light = light;
    m_transform = t;
    m_buf = Buffer3D(m_width, m_height, m_transform);
    findEyeRay();
    if (m_width > 0 && m_height > 0) {
        startThread();
    }
}

void Graph3D::findEyeRay() {
    bool invertible;
    Transform3D inv = m_transform.inverted(&invertible);
    Q_ASSERT(invertible);
    Vector3D pt1(0,0,1);
    Vector3D pt2(0,0,0);
    pt1 = inv * pt1;
    pt2 = inv * pt2;
    m_eyeray = (pt2 - pt1).normalized();
    qDebug() << m_eyeray << (inv * Vector3D(148, 149, -10) - inv * Vector3D(148, 149, 7)).normalized();
//    Q_ASSERT(qFuzzyCompare((inv * Vector3D(148, 149, -10) - inv * Vector3D(148, 149, 7)).normalized(), m_eyeray));
}

Sphere::Sphere(QObject* parent): Graph3D(parent) {
}

void Sphere::cancel() {
}

void Sphere::startThread() {
    for (float theta = 0; theta <= M_PI; theta += 0.01f) {
        for (float phi = 0; phi <= M_PI*2; phi += 0.01f) {
            Vector3D p(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));
            m_buf.drawTransformLitPoint(p, qRgba(0, 200, 0, 255), p, Vector3D(1, 1, 1));
        }
    }
    emit updated();
}

ImplicitGraph3D::ImplicitGraph3D(QObject* parent): Graph3D(parent), tv(Variable::Id("t", Variable::Id::Vector, te)) {
    constexpr int num = sizeof(te)/sizeof(te[0]);
    for (int i = 0; i < num; ++i) {
        te[i] = Number(i) / Number(num);
    }
}

void ImplicitGraph3D::cancel() {
    cancelled = true;
    future.waitForFinished();
    cancelled = false;
}

void ImplicitGraph3D::reset(std::unique_ptr<Equation> rel, const Variable& _x, const Variable& _y, const Variable& _z) {
    x = _x;
    y = _y;
    z = _z;
    func = Sub::create(std::move(rel->a), std::move(rel->b))->simplify();
    Expression::Subst s;
    EPtr ex = Variable::create(Variable::Id("v1x", Variable::Id::Constant, &v1[0])) + Variable::create(Variable::Id("dvx", Variable::Id::Constant, &dv[0])) * tv.ecopy();
    EPtr ey = Variable::create(Variable::Id("v1y", Variable::Id::Constant, &v1[1])) + Variable::create(Variable::Id("dvy", Variable::Id::Constant, &dv[1])) * tv.ecopy();
    EPtr ez = Variable::create(Variable::Id("v1z", Variable::Id::Constant, &v1[2])) + Variable::create(Variable::Id("dvz", Variable::Id::Constant, &dv[2])) * tv.ecopy();
    s.insert(std::make_pair(x, &*ex));
    s.insert(std::make_pair(y, &*ey));
    s.insert(std::make_pair(z, &*ez));
    rayfunc[0] = func->substitute(s)->simplify();
    rayfunc[1] = rayfunc[0]->derivative(tv)->simplify();
    d_rayfunc = (rayfunc[0]->ecopy() / rayfunc[1]->ecopy())->simplify();
    dx = func->derivative(x)->simplify();
    dy = func->derivative(y)->simplify();
    dz = func->derivative(z)->simplify();
}

void ImplicitGraph3D::startThread() {
    future = QtConcurrent::run(this, &ImplicitGraph3D::restart);
}

// [a, b)
template<typename T> QList<T> range(T a, T b) {
    QList<T> r;
    while (a < b) r.append(a++);
    return r;
}

void ImplicitGraph3D::restart() {
    bool invertible;
    const Transform3D inv = m_transform.inverted(&invertible);
    Q_ASSERT(invertible);
    QRgb c = m_color.rgba();
    for (int Y = 0; Y < m_height; ++Y) {
        UVector ox(VECTOR_ALLOC(m_width));
        UVector oy(VECTOR_ALLOC(m_width));
        UVector oz(VECTOR_ALLOC(m_width));
        std::unique_ptr<int[]> opt(new int[m_width]);
        std::size_t num = 0;
        for (int X = 0; X < m_width; ++X) {
            if (renderPoint(inv, Y, X, &ox[num], &oy[num], &oz[num])) {
                opt[num] = m_width*Y+X;
                ++num;
            }
            if (cancelled) return;
        }
        if (num) {
            x.id->type = Variable::Id::Vector;
            y.id->type = Variable::Id::Vector;
            z.id->type = Variable::Id::Vector;
            x.id->p = ox.get();
            y.id->p = oy.get();
            z.id->p = oz.get();
            VectorR vdx = dx->evaluateVector(num);
            VectorR vdy = dy->evaluateVector(num);
            VectorR vdz = dz->evaluateVector(num);
            for (std::size_t i = 0; i < num; ++i) {
                m_buf.drawTransformLitPoint(Vector3D(ox[i], oy[i], oz[i]), c, Vector3D(vdx[i], vdy[i], vdz[i]), m_light, opt[i]);
                if (cancelled) return;
            }
            VECTOR_FREE(vdx);
            VECTOR_FREE(vdy);
            VECTOR_FREE(vdz);
        }
        emit updated();
    }
}

template<int t> struct coord_tag {
    enum {
        index = t
    };
};
typedef coord_tag<0> X;
typedef coord_tag<1> Y;
typedef coord_tag<2> Z;
template<typename T> float get(const Vector3D& v) { return v.v.m[T::index]; }
template<typename T> void set(Vector3D& v, float n) { v.v.m[T::index] = n; }

struct Ray {
    Vector3D a, b;
    Ray() {}
    Ray(Vector3D _a, Vector3D _b) : a(_a), b(_b) { }
    Vector3D eval(float t) { return a + (b-a) * t; }
    template<typename T> float evalC(float t) const { return get<T>(a) + (get<T>(b) - get<T>(a)) * t; }
    template<typename T> EPtr evaluator(EPtr t) const {
        return Constant::create(get<T>(a)) + Constant::create(get<T>(b) - get<T>(a)) * std::move(t);
    }
};

template<typename T, typename U, typename V> bool intersect(const Ray& r, float t, float u1, float u2, float v1, float v2, float& q, Vector3D& pt) {
    float qn = t - get<T>(r.a), qd = get<T>(r.b) - get<T>(r.a);
    if (qd == 0) return false;
    q = qn / qd;
    float u = r.evalC<U>(q);
    if (u < u1 || u > u2) return false;
    float v = r.evalC<V>(q);
    if (v < v1 || v > v2) return false;
    set<T>(pt, t);
    set<U>(pt, u);
    set<V>(pt, v);
    return true;
}
template<typename T, typename U, typename V> bool intersect2(const Ray& r, const Vector3D& boxa, const Vector3D& boxb, float* q, Vector3D* pt, int& n) {
    (intersect<T, U, V>(r, get<T>(boxa), get<U>(boxa), get<U>(boxb), get<V>(boxa), get<V>(boxb), q[n], pt[n]) && ++n);
    if (n==2) {
        if (q[1]==q[0]) --n;
        else return true;
    }
    (intersect<T, U, V>(r, get<T>(boxb), get<U>(boxa), get<U>(boxb), get<V>(boxa), get<V>(boxb), q[n], pt[n]) && ++n);
    if (n==2) {
        if (q[1]==q[0]) --n;
        else return true;
    }
    return false;
}

inline bool rayAtPoint(const Transform3D& inv, const Vector3D& eyeray, float y, float x, const Vector3D& boxa, const Vector3D& boxb, Ray& ret) {
    Vector3D pp = inv * Vector3D(x, y, 1);
    Ray r(pp, pp + eyeray);
    int n = 0;
    Vector3D ends[2];
    float q[2];
    if (!(intersect2<X, Y, Z>(r, boxa, boxb, q, ends, n) || intersect2<Y, X, Z>(r, boxa, boxb, q, ends, n) || intersect2<Z, X, Y>(r, boxa, boxb, q, ends, n))) {
        // ray goes out of bounding box
        return false;
    }
    if (q[0] > q[1]) {
        // swap
        ret.a = ends[1];
        ret.b = ends[0];
    } else {
        ret.a = ends[0];
        ret.b = ends[1];
    }
    return true;
}

constexpr Number epsilon = 1e-5f;

inline bool zero(float n) { return !(n > epsilon || n < -epsilon); }

template<typename T> bool propernewton(const EPtr& f, const EPtr& g, const Variable& tv, Number& guess, T& guessChanged) {
    int iterations = 6;
    while (iterations--) {
        guess -= f->evaluate()/g->evaluate();
        guessChanged(guess);
    }
    Number v = f->evaluate();
    return gsl_finite(guess) && gsl_finite(v) && zero(v);
}

template<typename T> bool superbrute(const EPtr& f, const EPtr& g, const EPtr& h, const Variable& tv, Number* t, int size, Number& guess, T guessChanged) {
    UVector v(f->evaluateVector(size));
    UVector w(g->evaluateVector(size));
    int last = 0;
    for (int i = 0; i < size; ++i) {
        if (!gsl_finite(v[i])) {
            continue;
        }
#define DOGUESS \
        { \
            Number _t0 = t[0]; \
            t[0] = guess; \
            guessChanged(guess); \
            bool ok = propernewton(g, h, tv, t[0], guessChanged) && (qAbs(guess - t[0]) < 0.0625f) && t[0] >= 0 && t[0] <= 1; \
            guess = t[0]; \
            t[0] = _t0; \
            if (ok) { \
                return true; \
            } \
        }
        if (w[i] < epsilon && w[i] > -epsilon) {
            guess = v[i];
            DOGUESS
        }
        if (v[i]*v[last] < 0) {
            guess = (t[i]+t[last])/2;
            DOGUESS
        }
#undef DOGUESS
        last = i;
    }
    return false;
}

struct nullfunc { void operator()(...) { } };

bool ImplicitGraph3D::renderPoint(const Transform3D& inv, int py, int px, Vector ox, Vector oy, Vector oz) {
    Ray ray;
    if (!rayAtPoint(inv, m_eyeray, py, px, m_boxa, m_boxb, ray)) {
        return false;
    }
    v1[0] = ray.a.x();
    v1[1] = ray.a.y();
    v1[2] = ray.a.z();
    Vector3D _dv = ray.b - ray.a;
    dv[0] = _dv.x();
    dv[1] = _dv.y();
    dv[2] = _dv.z();
    Number guess = 0;
    if (!superbrute(d_rayfunc, rayfunc[0], rayfunc[1], tv, te, sizeof(te)/sizeof(te[0]), guess, nullfunc())) {
        return false;
    }
    Vector3D pt = ray.eval(guess);
    *ox = pt.x();
    *oy = pt.y();
    *oz = pt.z();
    return true;
}

QImage ImplicitGraph3D::diagnostics(const Transform3D& inv, int py, int px, Vector ox, Vector oy, Vector oz, QSize size) {
    Ray ray;
    if (!rayAtPoint(inv, m_eyeray, py, px, m_boxa, m_boxb, ray)) {
        return QImage();
    }
    QImage result(size, QImage::Format_RGB32);
    result.fill(Qt::white);
    v1[0] = ray.a.x();
    v1[1] = ray.a.y();
    v1[2] = ray.a.z();
    Vector3D _dv = ray.b - ray.a;
    dv[0] = _dv.x();
    dv[1] = _dv.y();
    dv[2] = _dv.z();
    int diagw = size.width(), diagh = size.height();
    constexpr int tesize = sizeof(te)/sizeof(te[0]);
    std::cerr << rayfunc[0]->toString() << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cerr << &v1[i] << ':' << v1[i] << std::endl;
    }
    for (int i = 0; i < 3; ++i) {
        std::cerr << &dv[i] << ':' << dv[i] << std::endl;
    }
    UVector f(rayfunc[0]->evaluateVector(tesize));
    UVector g(d_rayfunc->evaluateVector(tesize));
    float min = f[0], max = f[0];
    for (int i = 1; i < tesize; ++i) {
        min = qMin(min, f[i]);
        max = qMax(max, f[i]);
    }
    for (int i = 0; i < tesize; ++i) {
        result.setPixel(i*diagw/tesize, diagh-((-min)*diagh/(max-min)), qRgba(0, 0, 0, 255));
    }
    for (int i = 0; i < tesize; ++i) {
        result.setPixel(i*diagw/tesize, diagh-((f[i]-min)*diagh/(max-min)), qRgba(0, 0, 255, 255));
        result.setPixel(i*diagw/tesize, diagh-((g[i]-min)*diagh/(max-min)), qRgba(100, 100, 100, 255));
    }
    int wt = 0;
    auto fcn = std::function<void(float)>([&result,diagw,diagh,&wt](float g)->void{
        qDebug() << g;
        for (int yoy = wt; yoy < diagh; yoy += 8) {
            result.setPixel(g*diagw, yoy, qRgba(0, 255, 0, 255));
        }
        wt = ((wt+1)&7);
    });
    Number guess = 0;
    if (superbrute(d_rayfunc, rayfunc[0], rayfunc[1], tv, te, sizeof(te)/sizeof(te[0]), guess, fcn)) {
        fcn(guess);
        qDebug() << "yey" << guess;
    } else {
        qDebug() << "not okay" << guess;
    }
    return result;
}
