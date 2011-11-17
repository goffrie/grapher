#include "Graph.h"

#include <QtConcurrentMap>
#include <QtConcurrentRun>

#include <QPainter>

#include <QtGlobal>
#include <QDebug>
#include <iostream>
#include <functional>

#include <gsl/gsl_sys.h>
#include <gsl/gsl_nan.h>
#include <QTime>

constexpr int supersample = 2;

QImage downsample(QImage in) {
    return in.scaled(in.width() / supersample, in.height() / supersample, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
}

Graph::Graph(QObject* parent): QObject(parent) {
}

IteratingGraph::IteratingGraph(QObject* parent): Graph(parent), watcher(new QFutureWatcher<QImage>(this)) {
    connect(watcher, SIGNAL(finished()), this, SLOT(iterateAgain()));
}

void IteratingGraph::cancel() {
    cancelled = true;
    future.waitForFinished();
    cancelled = false;
}

void IteratingGraph::setupRestart(const QTransform& t, int _width, int _height) {
    cancel();
    width = _width * supersample;
    height = _height * supersample;
    transform = t;
    future = QtConcurrent::run(this, &IteratingGraph::restart);
    watcher->setFuture(future);
}

InequalityGraph::InequalityGraph(QObject* parent): Graph(parent) {
}

QImage InequalityGraph::img()
{
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

void InequalityGraph::setupRestart(const QTransform& t, int _width, int _height) {
    cancel();
    width = _width * supersample;
    height = _height * supersample;
    transform = t;
    future = QtConcurrent::run(this, &InequalityGraph::restart);
}

void InequalityGraph::restart() {
    QList<int> lines;
    for (int i = 0; i < height; ++i) lines.append(i);
    Number left, right, top, bottom;
    {
        QTransform ti = transform.inverted();
        QPointF tl = QPointF(0, 0) * ti;
        QPointF br = (QPointF(width-1, height-1) * 0.5) * ti;
        left = tl.x(); top = tl.y();
        right = br.x(); bottom = br.y();
    }
    Number xstep = (right - left) / width;
    Number ystep = (top - bottom) / height;
    QImage _img(width, height, QImage::Format_ARGB32_Premultiplied);
    std::function<void(int)> function([this, left, top, xstep, ystep, width, height, &_img](int y) -> void {
        if (this->cancelled) return;
        VectorR px = reinterpret_cast<VectorR>(alloca(sizeof(Number) * width));
        {
            Number _x = left;
            for (int x = 0; x < width; ++x) {
                px[x] = _x;
                _x += xstep;
            }
        }
        if (this->cancelled) return;
        Expression::Subst s;
        External X(px);
        Constant Y(top - y * ystep);
        s.insert(std::make_pair(this->x, &X));
        s.insert(std::make_pair(this->y, &Y));
        std::unique_ptr<Number[]> result(this->rel->substitute(s)->evaluateVector(width));
        if (this->cancelled) return;
        QRgb* __restrict__ p = reinterpret_cast<QRgb*>(_img.scanLine(y));
        for (int x = 0; x < width; ++x) {
            *p++ = qRgba(0, 0, 0, result[x] ? 255 : 0);
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
    cancelled = false;
}

ImplicitGraph::ImplicitGraph(QObject* parent) : IteratingGraph(parent) {
}

void ImplicitGraph::reset(const Equation& rel, const Variable& _x, const Variable& _y) {
    cancel();
    numPts = 0;
    std::cerr << rel.toString() << std::endl;
    eqn = Sub::create(rel.a->ecopy(), rel.b->ecopy());
    std::cerr << eqn->toString() << std::endl;
    x = _x;
    y = _y;
    _dx = eqn->derivative(x);
    _dy = eqn->derivative(y);
    std::cerr << _dx->toString() << '|' << _dy->toString() << std::endl;
    eqn = eqn->simplify();
    _dx = _dx->simplify();
    _dy = _dy->simplify();
    std::cerr << eqn->toString() << std::endl;
    std::cerr << _dx->toString() << '|' << _dy->toString() << std::endl;
    resubstitute();
}

QImage ImplicitGraph::restart() {
    QTransform ti = transform.inverted();
    m_px.reset(new Number[width * height]);
    m_py.reset(new Number[width * height]);
    resubstitute();
    numPts = 0;
    std::size_t size = 0;
    for (int y = 0; y < height; y += 3) {
        for (int x = 0; x < width; x += 3) {
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
    std::unique_ptr<Number[]> gx(dx->evaluateVector(size));
    if (cancelled) return QImage();
    std::unique_ptr<Number[]> gy(dy->evaluateVector(size));
    if (cancelled) return QImage();
    std::unique_ptr<Number[]> p (sub->evaluateVector(size));
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
    QImage _img(width, height, QImage::Format_ARGB32_Premultiplied);
    _img.fill(qRgba(0, 0, 0, 0));
    QPainter painter(&_img);
    painter.setPen(Qt::NoPen);
    painter.setBrush(Qt::black);
    painter.scale(supersample, supersample);
    for (std::size_t i = 0; i < numPts; ++i) {
        painter.fillRect(QRectF(QPointF(m_px[i], m_py[i]) * transform, QSizeF(0, 0)).adjusted(-0.5, -0.5, 0.5, 0.5), Qt::black);
    }
    painter.end();
    return _img;
}

void ImplicitGraph::resubstitute() {
    Expression::Subst s;
    External X(m_px.get()), Y(m_py.get());
    s.insert(std::make_pair(x, &X));
    s.insert(std::make_pair(y, &Y));
    dx = _dx->substitute(s);
    dy = _dy->substitute(s);
    sub = eqn->substitute(s);
}

QImage ImplicitGraph::iterate() {
    if (cancelled) return QImage();
    std::size_t size = numPts;
    std::unique_ptr<Number[]> p_gx (dx->evaluateVector(size));
    if (cancelled) return QImage();
    std::unique_ptr<Number[]> p_gy (dy->evaluateVector(size));
    if (cancelled) return QImage();
    std::unique_ptr<Number[]> p_p  (sub->evaluateVector(size));
    if (cancelled) return QImage();
    std::unique_ptr<Number[]> p_gx2(new Number[size]),
                              p_gy2(new Number[size]),
                              p_gxy(new Number[size]);
    VectorR gx  = p_gx.get(),
            gy  = p_gy.get(),
            p   = p_p.get(),
            gx2 = p_gx2.get(),
            gy2 = p_gy2.get(),
            gxy = p_gxy.get(),
            px  = m_px.get(),
            py  = m_py.get();
    // correcting factor = (gx, gy) * p / (gx * gx + gy * gy)
    // (gx2, gy2) = (gx^2, gy^2) in place
    for (std::size_t i = 0; i < size; ++i) {
        gx2[i] = gx[i] * gx[i];
    }
    if (cancelled) return QImage();
    for (std::size_t i = 0; i < size; ++i) {
        gy2[i] = gy[i] * gy[i];
    }
    if (cancelled) return QImage();
    // gxy = gx^2+gy^2
    for (std::size_t i = 0; i < size; ++i) {
        gxy[i] = gx2[i] + gy2[i];
    }
    if (cancelled) return QImage();
    // replace "p" with p/(gx^2 + gy^2)
    for (std::size_t i = 0; i < size; ++i) {
        p[i] = p[i] / gxy[i];
    }
    if (cancelled) return QImage();
    // replace (gx, gy) with (gx, gy) * p / (gx^2 + gy^2)
    for (std::size_t i = 0; i < size; ++i) {
        gx[i] = gx[i] * p[i];
    }
    if (cancelled) return QImage();
    for (std::size_t i = 0; i < size; ++i) {
        gy[i] = gy[i] * p[i];
    }
    if (cancelled) return QImage();
    // update points
    for (std::size_t i = 0; i < size; ++i) {
        px[i] = px[i] - gx[i];
    }
    if (cancelled) return QImage();
    for (std::size_t i = 0; i < size; ++i) {
        py[i] = py[i] - gy[i];
    }
    if (cancelled) return QImage();
    return downsample(draw());
}

ParametricGraph::ParametricGraph(QObject* parent): IteratingGraph(parent) {
}

void ParametricGraph::reset(std::unique_ptr<Expression> _x, std::unique_ptr<Expression> _y, const Variable& _t, Number _tMin, Number _tMax) {
    x = _x->simplify();
    y = _y->simplify();
    t = _t;
    tMin = _tMin;
    tMax = _tMax;
}

QImage ParametricGraph::restart() {
    const std::size_t num = 16;
    numPts = num;

    VectorR pt = new Number[num];
    EPtr sx, sy;
    {
        Expression::Subst s;
        External T(pt);
        s.insert(std::make_pair(t, &T));
        sx = x->substitute(s);
        sy = y->substitute(s);
    }
    Number step = (tMax - tMin) / (num - 1);
    Number n = tMin;
    for (std::size_t i = 0; i < num; ++i) {
        pt[i] = n;
        n += step;
    }
    VectorR vx, vy;
    {
        QFuture<Vector> fx = QtConcurrent::run(sx.get(), &Expression::evaluateVector, num);
        QFuture<Vector> fy = QtConcurrent::run(sy.get(), &Expression::evaluateVector, num);
        fx.waitForFinished();
        fy.waitForFinished();
        vx = fx.result();
        vy = fy.result();
    }
    _img = QImage(width, height, QImage::Format_ARGB32_Premultiplied);
    _img.fill(qRgba(0, 0, 0, 0));
    draw(vx, vy, num);
    m_vx.reset(vx);
    m_vy.reset(vy);
    m_pt.reset(pt);
    return downsample(_img);
}

void ParametricGraph::draw(Vector vx, Vector vy, size_t n) {
    QPainter painter(&_img);
    painter.setPen(Qt::NoPen);
    painter.setBrush(Qt::black);
    painter.scale(supersample, supersample);
    for (std::size_t i = 0; i < n; ++i) {
        painter.fillRect(QRectF(QPointF(vx[i], vy[i]) * transform, QSizeF(0, 0)).adjusted(-0.5, -0.5, 0.5, 0.5), Qt::black);
    }
    painter.end();
}

QImage ParametricGraph::iterate() {
    if (cancelled) return QImage();
    std::unique_ptr<Number[]> pt(new Number[numPts - 1]);
    std::unique_ptr<Number[]> nt(new Number[numPts * 2 - 1]);
    std::unique_ptr<Number[]> nx(new Number[numPts * 2 - 1]);
    std::unique_ptr<Number[]> ny(new Number[numPts * 2 - 1]);
    std::size_t numT = 0, numNT = 0;
    Number xscale = transform.m11(), yscale = transform.m22();
    bool kept = false;
    for (std::size_t i = 1; i < numPts; ++i) {
        Number dx = (m_vx[i] - m_vx[i-1]) * xscale, dy = (m_vy[i] - m_vy[i-1]) * yscale;
        if (gsl_finite(m_pt[i]) && gsl_finite(m_pt[i-1]) && m_pt[i] != m_pt[i-1] && dx*dx + dy*dy > 0.25) {
            // deepen
            if (!kept) {
                nt[numNT] = m_pt[i-1];
                nx[numNT] = m_vx[i-1];
                ny[numNT] = m_vy[i-1];
                ++numNT;
            }
            nt[numNT] = pt[numT] = (m_pt[i] + m_pt[i-1]) / 2;
            nx[numNT] = ny[numNT] = GSL_NAN;
            ++numNT;
            ++numT;
            
            nt[numNT] = m_pt[i];
            nx[numNT] = m_vx[i];
            ny[numNT] = m_vy[i];
            ++numNT;
            kept = true;
        } else {
            if (kept) {
                // fill in gaps with NaN to reduce wasted calculations
                nt[numNT] = GSL_NAN;
                nx[numNT] = ny[numNT] = GSL_NAN;
                ++numNT;
            }
            kept = false;
        }
    }
    if (cancelled) return QImage();
    qDebug() << "deepened by" << numT << "old" << numPts << "new" << numNT;
    if (numT == 0) {
        // done.
        m_pt.reset();
        m_vx.reset();
        m_vy.reset();
        numPts = 0;
        return downsample(_img);
    }
    std::unique_ptr<Number[]> vx, vy;
    {
        EPtr sx, sy;
        {
            Expression::Subst s;
            External T(pt.get());
            s.insert(std::make_pair(t, &T));
            sx = x->substitute(s);
            sy = y->substitute(s);
        }
        QFuture<Vector> fx = QtConcurrent::run(sx.get(), &Expression::evaluateVector, numT);
        QFuture<Vector> fy = QtConcurrent::run(sy.get(), &Expression::evaluateVector, numT);
        fx.waitForFinished();
        fy.waitForFinished();
        vx.reset(fx.result());
        vy.reset(fy.result());
    }
    if (cancelled) return QImage();
    QFuture<void> drawer = QtConcurrent::run(this, &ParametricGraph::draw, vx.get(), vy.get(), numT);
    for (std::size_t i = 0, j = 0; i < numT && j < numNT; ++j) {
        if (gsl_isnan(nt[j])) continue;
        if (nt[j] == pt[i]) {
            if (!gsl_isnan(ny[j])) {
                Q_ASSERT(!gsl_isnan(nx[j]));
                // problems.
                continue;
            }
            Q_ASSERT(gsl_isnan(nx[j]));
            if (!gsl_finite(vx[i]) || !gsl_finite(vy[i])) {
                continue;
            }
            nx[j] = vx[i];
            ny[j] = vy[i];
            ++i;
        }
    }
    numPts = numNT;
    m_vx = std::move(nx);
    m_vy = std::move(ny);
    m_pt = std::move(nt);
    drawer.waitForFinished();
    if (cancelled) return QImage();
    return downsample(_img);
}


void ParametricGraph::iterateAgain() {
    future.waitForFinished();
    if (future.isCanceled()) return;
    QImage newImg = future.result();
    if (newImg.isNull()) return;
    m_img = newImg;
    if (numPts != 0) {
        future = QtConcurrent::run(this, &ParametricGraph::iterate);
        watcher->setFuture(future);
    }
    emit updated();
}