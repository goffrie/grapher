#include "Graph2D.h"

#include <QtConcurrentMap>
#include <QtConcurrentRun>

#include <QPainter>

#include <QtGlobal>
#include <QDebug>
#include <iostream>
#include <functional>

#include <gsl/gsl_sys.h>
#include <gsl/gsl_nan.h>

#include <xmmintrin.h>
typedef __m128 v4sf;
#define VECTOR_LOOP(s) for (std::size_t i = 0; i < s; i += SSE_VECTOR_SIZE)
#define V(a) (*reinterpret_cast<v4sf*>(a+i))

inline QImage downsample(QImage in) {
    return in.scaled(in.width() / Graph2D::supersample, in.height() / Graph2D::supersample, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
}

Graph2D::Graph2D(QObject* parent): Graph(parent) {
}

void Graph2D::setupRestart(const QTransform& t, int _width, int _height) {
    cancel();
    m_width = _width * supersample;
    m_height = _height * supersample;
    transform = t;
    if (m_width > 0 && m_height > 0 && transform.isAffine()) {
        startThread();
    }
}

IteratingGraph::IteratingGraph(QObject* parent): Graph2D(parent), watcher(new QFutureWatcher<QImage>(this)) {
    connect(watcher, SIGNAL(finished()), this, SLOT(iterateAgain()));
}

void IteratingGraph::cancel() {
    cancelled = true;
    future.waitForFinished();
    future = QFuture<QImage>();
    cancelled = false;
}

void IteratingGraph::startThread() {
    future = QtConcurrent::run(this, &IteratingGraph::restart);
    watcher->setFuture(future);
}

InequalityGraph::InequalityGraph(QObject* parent): Graph2D(parent) {
}

QImage InequalityGraph::img() {
    img_mutex.lock();
    QImage r = m_img;
    img_mutex.unlock();
    return r;
}

void InequalityGraph::reset(std::unique_ptr<Inequality> _rel, const Variable& _x, const Variable& _y) {
    rel = _rel->simplify();
    x = _x;
    y = _y;
}

void InequalityGraph::startThread() {
    future = QtConcurrent::run(this, &InequalityGraph::restart);
}

void InequalityGraph::restart() {
    QList<int> lines;
    for (int i = 0; i < m_height; ++i) lines.append(i);
    Number left, right, top, bottom;
    {
        QTransform ti = transform.inverted();
        QPointF tl = QPointF(0, 0) * ti;
        QPointF br = (QPointF(m_width-1, m_height-1) * 0.5) * ti;
        left = tl.x(); top = tl.y();
        right = br.x(); bottom = br.y();
    }
    Number xstep = (right - left) / m_width;
    Number ystep = (top - bottom) / m_height;
    QImage _img(m_width, m_height, QImage::Format_ARGB32_Premultiplied);
    QRgb fill = qRgba(m_color.red(), m_color.green(), m_color.blue(), 255);
    std::function<void(int)> function([this, left, top, xstep, ystep, m_width, m_height, &_img, fill](int y) -> void {
        if (this->cancelled) return;
        VectorR px = reinterpret_cast<VectorR>(alloca(sizeof(Number) * m_width));
        {
            Number _x = left;
            for (int x = 0; x < m_width; ++x) {
                px[x] = _x;
                _x += xstep;
            }
        }
        if (this->cancelled) return;
        Expression::Subst s;
        Variable X(Variable::Id("X", Variable::Id::Vector, px));
        Constant Y(top - y * ystep);
        s.insert(std::make_pair(this->x, &X));
        s.insert(std::make_pair(this->y, &Y));
        UVector result(this->rel->substitute(s)->evaluateVector(m_width));
        if (this->cancelled) return;
        typedef QRgb* __restrict QRgbR;
        QRgbR p = reinterpret_cast<QRgb*>(_img.scanLine(y));
        for (int x = 0; x < m_width; ++x) {
            *p++ = result[x] ? fill : 0u;
        }
    });
    QtConcurrent::blockingMap(lines, function);
    if (cancelled) return;
    _img = downsample(_img);
    {
        QMutexLocker locker(&img_mutex);
        m_img = _img;
    }
    emit updated();
}

void InequalityGraph::cancel() {
    cancelled = true;
    future.waitForFinished();
    future = QFuture<void>();
    cancelled = false;
}

ImplicitGraph::ImplicitGraph(QObject* parent) : IteratingGraph(parent) {
}

void ImplicitGraph::reset(std::unique_ptr<Equation> rel, const Variable& _x, const Variable& _y) {
    cancel();
    numPts = 0;
    std::cerr << rel->toString() << std::endl;
    eqn = Sub::create(std::move(rel->a), std::move(rel->b))->simplify();
    std::cerr << eqn->toString() << std::endl;
    x = _x;
    y = _y;
    x.id->type = Variable::Id::Vector;
    y.id->type = Variable::Id::Vector;
    dx = eqn->derivative(x)->simplify();
    dy = eqn->derivative(y)->simplify();
    std::cerr << dx->toString() << '|' << dy->toString() << std::endl;
    resubstitute();
}

void ImplicitGraph::resubstitute() {
    x.id->p = m_px.get();
    y.id->p = m_py.get();
}

QImage ImplicitGraph::restart() {
    QTransform ti = transform.inverted();
    m_px.reset(VECTOR_ALLOC(m_width * m_height));
    m_py.reset(VECTOR_ALLOC(m_width * m_height));
    resubstitute();
    numPts = 0;
    std::size_t size = 0;
    for (int y = 0; y < m_height; y += 3) {
        for (int x = 0; x < m_width; x += 3) {
            {
            QPointF pt = QPointF(x, y) * ti;
            m_px[size] = pt.x();
            m_py[size] = pt.y();
            ++size;
            }
            {
            QPointF pt = QPointF(x + 2, y + 1) * ti;
            m_px[size] = pt.x();
            m_py[size] = pt.y();
            ++size;
            }
            {
            QPointF pt = QPointF(x + 1, y + 2) * ti;
            m_px[size] = pt.x();
            m_py[size] = pt.y();
            ++size;
            }
        }
        if (cancelled) return QImage();
    }
    UVector gx(dx->evaluateVector(size));
    if (cancelled) return QImage();
    UVector gy(dy->evaluateVector(size));
    if (cancelled) return QImage();
    UVector p (eqn->evaluateVector(size));
    if (cancelled) return QImage();
    for (std::size_t i = 0; i < size; ++i) {
        QPointF pt(m_px[i], m_py[i]);
        QPointF correction = QPointF(gx[i], gy[i]) * (p[i] / (gx[i] * gx[i] + gy[i] * gy[i]));
        correction = (pt * transform) - ((pt - correction) * transform);
        qreal len = correction.x() * correction.x() + correction.y() * correction.y();
        if (len < 200) {
            m_px[numPts] = pt.x();
            m_py[numPts] = pt.y();
            ++numPts;
        }
    }
    if (cancelled) return QImage();
    return downsample(draw());
}

void ImplicitGraph::iterateAgain() {
    future.waitForFinished();
    if (future.isCanceled()) return;
    QImage newImg = future.result();
    if (newImg.isNull()) return;
    m_img = newImg;
    future = QtConcurrent::run(this, &ImplicitGraph::iterate);
    watcher->setFuture(future);
    emit updated();
}

QImage ImplicitGraph::draw() {
    QImage _img(m_width, m_height, QImage::Format_ARGB32_Premultiplied);
    _img.fill(qRgba(0, 0, 0, 0));
    QPainter painter(&_img);
    if (!painter.isActive()) return QImage();
    painter.scale(supersample, supersample);
    for (std::size_t i = 0; i < numPts; ++i) {
        painter.fillRect(QRectF(QPointF(m_px[i], m_py[i]) * transform, QSizeF(0, 0)).adjusted(-0.5, -0.5, 0.5, 0.5), m_color);
    }
    painter.end();
    return _img;
}

QImage ImplicitGraph::iterate() {
    if (cancelled) return QImage();
    std::size_t size = numPts;
    UVector p_gx (dx->evaluateVector(size));
    if (cancelled) return QImage();
    UVector p_gy (dy->evaluateVector(size));
    if (cancelled) return QImage();
    UVector p_p  (eqn->evaluateVector(size));
    if (cancelled) return QImage();
    const VectorR
        gx = p_gx.get(),
        gy = p_gy.get(),
        p  = p_p.get(),
        px = m_px.get(),
        py = m_py.get();
    // correcting factor = (gx, gy) * p / (gx * gx + gy * gy)
    // (gx2, gy2) = (gx^2, gy^2) in place
    VECTOR_LOOP(size) {
        // replace "p" with p/(gx^2 + gy^2)
        v4sf _gx = V(gx), _gy = V(gy);
        v4sf _p = V(p) / (_gx * _gx + _gy * _gy);
        // delta-P = (gx, gy) * p / (gx^2 + gy^2)
        // update points
        V(px) -= _gx * _p;
        V(py) -= _gy * _p;
    }
    if (cancelled) return QImage();
    return downsample(draw());
}

ParametricGraph::ParametricGraph(QObject* parent): IteratingGraph(parent), pts(VECTOR_ALLOC(numPts)) {
}

void ParametricGraph::reset(std::unique_ptr<Expression> _x, std::unique_ptr<Expression> _y, const Variable& _t, Number _tMin, Number _tMax) {
    x = _x->simplify();
    y = _y->simplify();
    t = _t;
    t.id->type = Variable::Id::Vector;
    t.id->p = pts.get();
    tMin = _tMin;
    tMax = _tMax;
}

QImage ParametricGraph::restart() {
    if (tMin == 0 && tMax == 0) {
        Variable* vx = dynamic_cast<Variable*>(x.get());
        Variable* vy = dynamic_cast<Variable*>(y.get());
        QTransform ti = transform.inverted();
        if (vx && *vx == t) {
            tMin = (QPointF(0, 0) * ti).x();
            tMax = (QPointF(m_width/supersample, 0) * ti).x();
        } else if (vy && *vy == t) {
            tMax = (QPointF(0, 0) * ti).y();
            tMin = (QPointF(0, m_height/supersample) * ti).y();
        }
        Q_ASSERT(tMax > tMin);
    }
    distribution = std::uniform_real_distribution<Number>(tMin, tMax);

    _img = QImage(m_width, m_height, QImage::Format_ARGB32_Premultiplied);
    _img.fill(qRgba(0, 0, 0, 0));
    return iterate();
}

void ParametricGraph::draw(Vector vx, Vector vy) {
    QPainter painter(&_img);
    painter.scale(supersample, supersample);
    for (std::size_t i = 0; i < numPts; ++i) {
        painter.fillRect(QRectF(QPointF(vx[i], vy[i]) * transform, QSizeF(0, 0)).adjusted(-0.5, -0.5, 0.5, 0.5), m_color);
    }
    painter.end();
}

QImage ParametricGraph::iterate() {
    if (cancelled) return QImage();
    for (int i = 0; i < numPts; ++i) pts[i] = distribution(engine);
    if (cancelled) return QImage();
    QFuture<Vector> fy = QtConcurrent::run(y.get(), &Expression::evaluateVector, (std::size_t)numPts);
    UVector vx(x->evaluateVector(numPts));
    fy.waitForFinished();
    UVector vy(fy.result());
    if (cancelled) return QImage();
    draw(vx.get(), vy.get());
    if (cancelled) return QImage();
    return downsample(_img);
}

void ParametricGraph::iterateAgain() {
    future.waitForFinished();
    if (future.isCanceled()) return;
    QImage newImg = future.result();
    if (newImg.isNull()) return;
    m_img = newImg;
    future = QtConcurrent::run(this, &ParametricGraph::iterate);
    watcher->setFuture(future);
    emit updated();
}
