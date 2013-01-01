#include "ImplicitGraph3D.h"

#include "global.h"

#include <QPainter>

#include <iomanip>
#include <iostream>
#include <functional>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <complex>

#ifdef _WIN32
#include <malloc.h>
#else
#include <alloca.h>
#endif

#include <gsl/gsl_sys.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_math.h>


ImplicitGraph3D::ImplicitGraph3D(QObject* parent): Graph3D(parent), tv(Variable::Id("t", Variable::Id::Vector, te->begin())),
v1x(Variable::Id("v1x", Variable::Id::Constant, &v1[0])),
v1y(Variable::Id("v1y", Variable::Id::Constant, &v1[1])),
v1z(Variable::Id("v1z", Variable::Id::Constant, &v1[2])),
dvx(Variable::Id("dvx", Variable::Id::Constant, &dv[0])),
dvy(Variable::Id("dvy", Variable::Id::Constant, &dv[1])),
dvz(Variable::Id("dvz", Variable::Id::Constant, &dv[2]))
{
    auto& _te = *te;
    for (int i = 0; i < _te.size(); ++i) {
        _te[i] = Number(i) / _te.size();
    }
}

ImplicitGraph3D::~ImplicitGraph3D() {
}

// HACK: make this not a global variable, etc.
static MathContext* ctx = nullptr;

void ImplicitGraph3D::reset(std::unique_ptr<Equation> rel, const Variable& _x, const Variable& _y, const Variable& _z) {
    if (!ctx) ctx = new MathContext(MathContext::defaultContext());
    x = _x;
    y = _y;
    z = _z;
    func = Sub::create(std::move(rel->a), std::move(rel->b))->simplify();
    Expression::Subst s;
    EPtr ex = v1x.ecopy() + dvx.ecopy() * tv.ecopy();
    EPtr ey = v1y.ecopy() + dvy.ecopy() * tv.ecopy();
    EPtr ez = v1z.ecopy() + dvz.ecopy() * tv.ecopy();
    s.insert(std::make_pair(x, &*ex));
    s.insert(std::make_pair(y, &*ey));
    s.insert(std::make_pair(z, &*ez));
    rayfunc[0] = func->substitute(s)->simplify();
    if (rayfunc[0]->polynomial(tv)) {
        polyrayfunc = rayfunc[0]->expand()->facsum(tv);
        Polynomial* poly = polyrayfunc.get();
        int d = poly->degree();
        polyrayfunc_e.resize(d+1);
        for (int i = d; i >= 0; --i) {
            Q_ASSERT(poly);
            polyrayfunc_e[i] = poly->right->evaluator(*ctx);
            poly = poly->left.get();
        }
        Q_ASSERT(!poly);
    }
    rayfunc[1] = rayfunc[0]->derivative(tv)->simplify();
    rayfunc[2] = rayfunc[1]->derivative(tv)->simplify();
    for (int i = 0; i < 3; ++i) {
        rayfunc_e[i] = rayfunc[i]->evaluator(*ctx);
    }
    //rayfunc_v = rayfunc[0]->evaluatorVector();
    d_rayfunc = (rayfunc[0]->ecopy() / rayfunc[1]->ecopy())->simplify();
    //d_rayfunc_v = d_rayfunc->evaluatorVector();
    dx = func->derivative(x)->simplify();
    dy = func->derivative(y)->simplify();
    dz = func->derivative(z)->simplify();
}

// [a, b)
template<typename T> QList<T> range(T a, T b) {
    QList<T> r;
    while (a < b) r.append(a++);
    return r;
}

void ImplicitGraph3D::compute() {
    _mm_empty();
    const Transform3D inv = transform().inverted();
    Buffer3D buffer(width(), height(), transform());
    buffer.setColor(color().rgba());
    buffer.setLight(light());

    std::unique_ptr<int[]> opt(new int[width()]);
    UVector ox(VECTOR_ALLOC(width()));
    UVector oy(VECTOR_ALLOC(width()));
    UVector oz(VECTOR_ALLOC(width()));
    x.id->type = Variable::Id::Vector;
    y.id->type = Variable::Id::Vector;
    z.id->type = Variable::Id::Vector;
    x.id->p = ox.get();
    y.id->p = oy.get();
    z.id->p = oz.get();
    for (int Y = 0; Y < height(); ++Y) {
        std::size_t num = 0;
        for (int X = 0; X < width(); ++X) {
            if (renderPoint(inv, X, Y, &ox[num], &oy[num], &oz[num])) {
                opt[num] = width()*Y+X;
                ++num;
            }
            if (cancelled()) return;
        }
        if (num) {
            // compute normal vector
            UVector vdx(dx->evaluateVector(num));
            UVector vdy(dy->evaluateVector(num));
            UVector vdz(dz->evaluateVector(num));
            // draw points
            _mm_empty();
            for (std::size_t i = 0; i < num; ++i) {
                buffer.drawTransformLitPoint(Vector3D<float>(ox[i], oy[i], oz[i]), Vector3D<float>(vdx[i], vdy[i], vdz[i]), opt[i]);
                if (cancelled()) return;
            }
        }
        emit updated(new Buffer3D(buffer.copy()));
    }
}

BOOST_CONSTEXPR_OR_CONST int X = 0, Y = 1, Z = 2;

struct Ray {
    Vector3D<float> a, b;
    Ray() {}
    Ray(Vector3D<float> _a, Vector3D<float> _b) : a(_a), b(_b) { }
    Vector3D<float> eval(float t) { return a + (b-a) * t; }
    template<int T> float evalC(float t) const { return a.get<T>() + (b.get<T>() - a.get<T>()) * t; }
    template<int T> EPtr evaluator(EPtr t) const {
        return Constant::create(a.get<T>()) + Constant::create(b.get<T>() - a.get<T>()) * std::move(t);
    }
};

template<int T, int U, int V> bool intersect(const Ray& r, float t, float u1, float u2, float v1, float v2, float& q, Vector3D<float>& pt) {
    float qn = t - r.a.get<T>(), qd = r.b.get<T>() - r.a.get<T>();
    if (qd == 0) return false;
    q = qn / qd;
    float u = r.evalC<U>(q);
    if (u < u1 || u > u2) return false;
    float v = r.evalC<V>(q);
    if (v < v1 || v > v2) return false;
    pt.set<T>(t);
    pt.set<U>(u);
    pt.set<V>(v);
    return true;
}
template<int T, int U, int V> bool intersect2(const Ray& r, const Vector3D<float>& boxa, const Vector3D<float>& boxb, float* q, Vector3D<float>* pt, int& n) {
    (intersect<T, U, V>(r, boxa.get<T>(), boxa.get<U>(), boxb.get<U>(), boxa.get<V>(), boxb.get<V>(), q[n], pt[n]) && ++n);
    if (n==2) {
        if (q[1]==q[0]) --n;
        else return true;
    }
    (intersect<T, U, V>(r, boxb.get<T>(), boxa.get<U>(), boxb.get<U>(), boxa.get<V>(), boxb.get<V>(), q[n], pt[n]) && ++n);
    if (n==2) {
        if (q[1]==q[0]) --n;
        else return true;
    }
    return false;
}

inline bool rayAtPoint(const Transform3D& inv, const Vector3D<float>& eyeray, float y, float x, const Vector3D<float>& boxa, const Vector3D<float>& boxb, Ray& ret) {
    Vector3D<float> pp = inv * Vector3D<float>(x, y, 1);
    Ray r(pp, pp + eyeray);
    int n = 0;
    Vector3D<float> ends[2];
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

BOOST_CONSTEXPR_OR_CONST Number epsilon = 1.f / (1<<8);
BOOST_CONSTEXPR_OR_CONST Number bigepsilon = 1.f / (1<<4);

template<typename F> bool zero(F n) { return !(n > epsilon || n < -epsilon); }
inline bool zero_wrt(float n, float d) {
    return zero(n / (std::fabs(d)+1));
}

template<typename T> bool newton(const WEvalFunc& f, const WEvalFunc& g, const Variable& tv, Number& guess, T& guessChanged) {
    int iterations = 10;
    float d;
    while (iterations--) {
        guess -= (d = f())/g();
        guessChanged(guess);
    }
    Number v = f();
    return gsl_finite(guess) && gsl_finite(v) && zero_wrt(v, d);
}

template<bool diag, typename T> bool halley(const WEvalFunc& f, const WEvalFunc& g, const WEvalFunc& h, const Variable& tv, Number& guess, T& guessChanged) {
    int iterations = 5;
    float G;
    while (iterations--) {
        Number F = f(), H = h();
        G = g();
        guess -= F * G / (G * G - 0.5f * F * H);
        guessChanged(guess);
    }
    Number v = f();
    if (gsl_finite(guess)) {
        if (gsl_finite(v)) {
            if (zero_wrt(v, G)) {
                return true;
            } else if (diag) {
                float e = epsilon * (std::fabs(G)+1);
                qDebug() << "Rejected" << guess << "for f(t) =" << v << "being outside of [" << -e << "," << e << "]";
            }
        } else if (diag) {
            qDebug() << "Rejected" << guess << "for f(t) =" << v << "being non-finite";
        }
    } else if (diag) {
        qDebug() << "Rejected non-finite guess" << guess;
    }
    return false;
}

template<bool diag = false, typename T> bool superbrute(const Expression& f, const Expression& gv, const WEvalFunc& g, const WEvalFunc& h, const WEvalFunc& j, const Variable& tv, Number* t, int size, Number& guess, T guessChanged) {
    UVector v(f.evaluateVector(size)), w(gv.evaluateVector(size));
    int last = 0;
    const float maxdelta = 10.f / size;
    for (int i = 0; i < size; ++i) {
        if (!gsl_finite(v[i])) {
            continue;
        }
        bool go = true;
        if (v[i]*v[last] < 0) {
            guess = (t[i]+t[last])/2;
        } else if (w[i] < bigepsilon && w[i] > -bigepsilon) {
            guess = t[i];
        } else {
            go = false;
        }
        if (go) {
            Number _t0 = t[0];
            t[0] = guess;
            guessChanged(guess);
            bool ok = false;
            if (halley<diag>(g, h, j, tv, t[0], guessChanged)) {
                if (qAbs(guess - t[0]) < maxdelta) {
                    if (t[0] >= -epsilon && t[0] <= 1+epsilon) {
                        ok = true;
                    } else if (diag) {
                        qDebug() << "Rejected" << t[0] << "for being outside of [0, 1]";
                    }
                } else if (diag) {
                    qDebug() << "Rejected" << t[0] << "for being too far from initial guess" << guess;
                }
            }
            guess = t[0];
            t[0] = _t0;
            if (ok) {
                return true;
            }
        }
        last = i;
    }
    return false;
}

struct Poly {
    int degree;
    double* coeff;
    Poly(const std::vector<WEvalFunc>& poly_coeff) {
        degree = poly_coeff.size()-1;
        coeff = new double[degree];
        double factor = 1. / poly_coeff[0]();
        for (int i = 0; i < degree; ++i) {
            coeff[i] = factor * poly_coeff[i+1]();
        }
    }
    ~Poly() {
        delete[] coeff;
    }
    template<typename T> std::complex<T> evaluate(std::complex<T> x) const {
        if (degree == 0) return x;
        std::complex<T> r = x + coeff[0];
        for (int i = 1; i < degree; ++i) {
            r *= x;
            r += coeff[i];
        }
        return r;
    }
};
std::ostream& operator<<(std::ostream& s, const Poly& p) {
    s << "x^" << p.degree;
    for (int i = 0; i < p.degree; ++i) {
        s << std::setprecision(14) << " + " << p.coeff[i];
        if (p.degree-i-1 > 0) s << "*x";
        if (p.degree-i-1 > 1) s << '^' << (p.degree-i-1);
    }
    return s;
}

template<bool diag = false, typename T>
void durandkerner(const Poly& poly, std::complex<double>* roots, T rootChanged) {
    if (diag) std::cerr << poly << std::endl;
    int degree = poly.degree;
    if (degree == 0) {
    	return;
	}
    if (degree == 1) {
        // x + a = 0 => x = -a
        roots[0] = -poly.coeff[0];
        rootChanged(roots[0].real());
        return;
    }
    if (degree == 2) {
        // x^2 + ax + b = 0
        // x = (-a +/- sqrt(a^2 - 4b)) / 2
        std::complex<double> disc = std::sqrt(poly.coeff[0] * poly.coeff[0] - 4.0 * poly.coeff[1]);
        roots[0] = 0.5 * (-poly.coeff[0] + disc);
        rootChanged(roots[0].real());
        roots[1] = 0.5 * (-poly.coeff[0] - disc);
        rootChanged(roots[1].real());
        return;
    }
    BOOST_CONSTEXPR_OR_CONST std::complex<double> startbase(0.6, 0.94);
    roots[0] = 1.;
    for (int i = 1; i < degree; ++i) {
        roots[i] = roots[i-1] * startbase;
    }
    int maxiterations = 30;
    bool repeat = true;
    while (repeat && maxiterations--) {
        repeat = false;
        for (int i = 0; i < degree; ++i) {
            std::complex<double> denom;
            bool first = true;
            for (int j = 0; j < degree; ++j) { if (i == j) continue;
                if (first) {
                    denom = roots[i] - roots[j];
                    first = false;
                } else {
                    denom *= roots[i] - roots[j];
                }
            }
            const std::complex<double> adjust = poly.evaluate(roots[i]) / denom;
            if (adjust.real() * adjust.real() + adjust.imag() * adjust.imag() > 0.0000001) repeat = true;
            roots[i] -= adjust;
            if (diag) std::cerr << roots[i] << ' ';
            rootChanged(roots[i].real());
        }
        if (diag) std::cerr << std::endl;
    }
}

template<bool diag = false, typename T>
bool polyroot(const std::vector<WEvalFunc>& poly_coeff, float& root, T rootChanged) {
    Poly _poly(poly_coeff);
    const int degree = _poly.degree;
    std::complex<double>* croots = reinterpret_cast<std::complex<double>*>(alloca(sizeof(std::complex<double>)*degree));
    durandkerner<diag>(_poly, croots, rootChanged);
    bool ok = false;
    for (int i = 0; i < degree; ++i) {
        if (zero(croots[i].imag())) {
            float r = croots[i].real();
            if (gsl_finite(r) && r >= -epsilon && r <= 1+epsilon && (!ok || r < root)) {
                root = r;
                ok = true;
            }
        }
    }
    return ok;
}

struct nullfunc { void operator()(...) { } };

bool ImplicitGraph3D::renderPoint(const Transform3D& inv, int px, int py, Vector ox, Vector oy, Vector oz) {
    Ray ray;
    if (!rayAtPoint(inv, eyeray(), py, px, boxa(), boxb(), ray)) {
        return false;
    }
    v1[0] = ray.a.x();
    v1[1] = ray.a.y();
    v1[2] = ray.a.z();
    Vector3D<float> _dv = ray.b - ray.a;
    dv[0] = _dv.x();
    dv[1] = _dv.y();
    dv[2] = _dv.z();
    Number guess = 0;
    if (polyrayfunc.get()) {
        if (!polyroot(polyrayfunc_e, guess, nullfunc())) {
            return false;
        }
    } else {
        if (!superbrute(*d_rayfunc, *rayfunc[0], rayfunc_e[0], rayfunc_e[1], rayfunc_e[2], tv, te->begin(), te->size(), guess, nullfunc())) {
            return false;
        }
    }
    Vector3D<float> pt = ray.eval(guess);
    *ox = pt.x();
    *oy = pt.y();
    *oz = pt.z();
    return true;
}

inline void analyze(Vector data, int size, float& mean, float& stdev) {
    mean = 0;
    float factor = 1.f / size;
    for (int i = 0; i < size; ++i) {
        mean += data[i] * factor;
    }
    stdev = 0;
    for (int i = 0; i < size; ++i) {
        float resid = data[i] - mean;
        stdev += resid * resid * factor;
    }
    stdev = std::sqrt(stdev);
}
inline void quantile(Vector data, int size, int which, float& q1, float& q3) {
    Vector v = VECTOR_ALLOC(size);
    std::memcpy(v, data, sizeof(Number)*size);
    std::nth_element(v, v+size/which, v+size);
    std::nth_element(v, v+size*(which-1)/which, v+size);
    q1 = v[size/which];
    q3 = v[size*(which-1)/which];
    VECTOR_FREE(v);
}

QPixmap ImplicitGraph3D::diagnostics(const Transform3D& inv, int px, int py, QSize size) {
    QPixmap result(size);
    result.fill(Qt::white);
    Ray ray;
    if (!rayAtPoint(inv, eyeray(), py, px, boxa(), boxb(), ray)) {
        return result;
    }
    QPainter painter(&result);
    painter.fillRect(10, 10, 10, 10, Qt::blue);
    painter.scale(1, -1);
    painter.translate(1, -size.height()+1);
    painter.fillRect(10, 10, 10, 10, Qt::black);
    v1[0] = ray.a.x();
    v1[1] = ray.a.y();
    v1[2] = ray.a.z();
    Vector3D<float> _dv = ray.b - ray.a;
    dv[0] = _dv.x();
    dv[1] = _dv.y();
    dv[2] = _dv.z();
    {
        int diagw = size.width() - 2, diagh = size.height() - 2;
        painter.scale(diagw, diagh);
        painter.scale(1, 0.5);
        painter.translate(0, 1);
    }
    BOOST_CONSTEXPR_OR_CONST int tesize = std::tuple_size<std::remove_reference<decltype(*te)>::type>::value;
    Expression::Subst s;
    Constant _v1x(v1[0]);
    Constant _v1y(v1[1]);
    Constant _v1z(v1[2]);
    Constant _dvx(dv[0]);
    Constant _dvy(dv[1]);
    Constant _dvz(dv[2]);
    s.insert(std::make_pair(v1x, &_v1x));
    s.insert(std::make_pair(v1y, &_v1y));
    s.insert(std::make_pair(v1z, &_v1z));
    s.insert(std::make_pair(dvx, &_dvx));
    s.insert(std::make_pair(dvy, &_dvy));
    s.insert(std::make_pair(dvz, &_dvz));
    {
        EPtr p = rayfunc[0]->substitute(s);
        std::cerr << p->toString() << std::endl << std::endl;
        p = p->expand();
        std::cerr << p->toString() << std::endl << std::endl;
        p = p->facsum(tv);
        std::cerr << p->toString() << std::endl << std::endl;
        std::cerr << '!' << p->simplify()->toString() << std::endl << std::endl;
    }
    if (polyrayfunc.get()) {
        std::cerr << polyrayfunc->toString() << std::endl << std::endl;
        //std::cerr << polyrayfunc->substitute(s)->toString() << std::endl;
        std::cerr << polyrayfunc->substitute(s)->facsum(tv)->toString() << std::endl;
        std::cerr << '?' << polyrayfunc->substitute(s)->facsum(tv)->simplify()->toString() << std::endl;
        for (uz i = 0; i < polyrayfunc_e.size(); ++i) {
            std::cerr << polyrayfunc_e[i]() << ' ';
        }
        std::cerr << std::endl;
    }
    UVector f(rayfunc[0]->evaluateVector(tesize));
    UVector g(d_rayfunc->evaluateVector(tesize));
    painter.setPen(Qt::black);
    painter.drawLine(0, 0, 1, 0);
    {
        painter.save();
/*        float mean, stdev;
        analyze(f.get(), tesize, mean, stdev);
        float min = mean - stdev, max = mean + stdev;
        painter.scale(1, 1.f/qMax(qAbs(max), qAbs(min)));*/
        float q1, q3;
        quantile(f.get(), tesize, 12, q1, q3);
        painter.scale(1, 0.5f/qMin(qAbs(q1), qAbs(q3)));
        painter.setPen(QColor(0, 0, 255));
        QPointF lastPoint(0, GSL_NAN);
        for (int i = 0; i < tesize; ++i) {
            QPointF newPoint((*te)[i], f[i]);
            if (gsl_finite(lastPoint.y())) {
                painter.drawLine(lastPoint, newPoint);
            } else {
                painter.drawPoint(newPoint);
            }
            lastPoint = newPoint;
        }
//        result.setPixel(i*diagw/tesize, diagh-((f[i]-min)*diagh/(max-min)));
        painter.restore();
    }
    if (polyrayfunc) {
    	UVector r(polyrayfunc->evaluateVector(tesize));
        painter.save();
/*        float mean, stdev;
        analyze(f.get(), tesize, mean, stdev);
        float min = mean - stdev, max = mean + stdev;
        painter.scale(1, 1.f/qMax(qAbs(max), qAbs(min)));*/
        float q1, q3;
        quantile(f.get(), tesize, 12, q1, q3);
        painter.scale(1, 0.5f/qMin(qAbs(q1), qAbs(q3)));
        painter.setPen(QColor(200, 0, 200));
        QPointF lastPoint(0, GSL_NAN);
        for (int i = 0; i < tesize; ++i) {
            QPointF newPoint((*te)[i], r[i]);
            if (gsl_finite(lastPoint.y())) {
                painter.drawLine(lastPoint, newPoint);
            } else {
                painter.drawPoint(newPoint);
            }
            lastPoint = newPoint;
        }
//        result.setPixel(i*diagw/tesize, diagh-((f[i]-min)*diagh/(max-min)));
        painter.restore();
    }
    {
        painter.save();
/*        float mean, stdev;
        analyze(g.get(), tesize, mean, stdev);
        float min = mean - stdev, max = mean + stdev;
        painter.scale(1, 1.f/qMax(qAbs(max), qAbs(min)));*/
        float q1, q3;
        quantile(g.get(), tesize, 12, q1, q3);
        painter.scale(1, 0.5f/qMin(qAbs(q1), qAbs(q3)));
        painter.setPen(QColor(100, 100, 100));
        QPointF lastPoint(0, GSL_NAN);
        for (int i = 0; i < tesize; ++i) {
            QPointF newPoint((*te)[i], g[i]);
            if (gsl_finite(lastPoint.y())) {
                painter.drawLine(lastPoint, newPoint);
            } else {
                painter.drawPoint(newPoint);
            }
            lastPoint = newPoint;
        }
//        result.setPixel(i*diagw/tesize, diagh-((g[i]-min)*diagh/(max-min)));
        painter.restore();
    }
    int wt = 0;
    painter.setPen(QColor(0, 255, 0, 100));
    auto fcn = std::function<void(float)>([&painter,&wt](float g)->void{
        qDebug() << g;
        painter.drawLine(QPointF(g, -1), QPointF(g, 1));
/*        for (int yoy = wt; yoy < diagh; yoy += 8) {
            result.setPixel(g*diagw, yoy, qRgba(0, 255, 0, 255));
        }*/
        wt = ((wt+1)&7);
    });
    Number guess = 0;
    if (polyrayfunc.get()) {
        if (polyroot<true>(polyrayfunc_e, guess, fcn)) {
            qDebug() << "yey" << guess;
        } else {
            qDebug() << "not okay" << guess;
        }
    } else {
        if (superbrute<true>(*d_rayfunc, *rayfunc[0], rayfunc_e[0], rayfunc_e[1], rayfunc_e[2], tv, te->begin(), te->size(), guess, fcn)) {
            fcn(guess);
            qDebug() << "yey" << guess;
        } else {
            qDebug() << "not okay" << guess;
        }
    }
    painter.end();
    return result;
}
