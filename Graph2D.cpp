#include "Graph2D.h"

#include <QtConcurrentRun>
#include <QPainter>
#include <QtGlobal>
#include <QDebug>

#include <gsl/gsl_sys.h>
#include <gsl/gsl_nan.h>

#include "global.h"

#include "Expression.h"

inline QImage downsample(QImage in) {
    return in.scaled(in.width() / Graph2D::supersample, in.height() / Graph2D::supersample, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
}

Graph2D::Graph2D(QObject* parent): Graph(parent) {
}

void Graph2D::setupRestart(const QTransform& transform, QSize size) {
    stop();
    m_size = size * supersample;
    m_transform = transform;
    if (m_size.width() > 0 && m_size.height() > 0 && m_transform.isAffine()) {
        restart();
    }
}

InequalityGraph::InequalityGraph(QObject* parent): Graph2D(parent) {
}

void InequalityGraph::reset(std::unique_ptr<Inequality> rel, const Variable& _x, const Variable& _y) {
    m_rel = rel->simplify();
    x = _x;
    y = _y;
}

void InequalityGraph::compute() {
    const int width = size().width(), height = size().height();
    Number left, right, top, bottom;
    {
        QTransform inverted = transform().inverted();
        QPointF topLeft = QPointF(0, 0) * inverted;
        QPointF bottomRight = (QPointF(width - 1, height - 1) * 0.5) * inverted;
        left = topLeft.x(); top = topLeft.y();
        right = bottomRight.x(); bottom = bottomRight.y();
    }
    Number xstep = (right - left) / size().width();
    Number ystep = (top - bottom) / size().height();
    QImage img(width, height, QImage::Format_ARGB32_Premultiplied);
    QRgb fill = color().rgb();
    UVector px(VECTOR_ALLOC(width));
    {
        Number _x = left;
        for (int x = 0; x < width; ++x) {
            px[x] = _x;
            _x += xstep;
        }
    }
    float lineY;
    x.id->type = Variable::Id::Vector;
    x.id->p = px.get();
    y.id->type = Variable::Id::Constant;
    y.id->p = &lineY;
    for (int line = 0; line < height; ++line) {
        if (cancelled()) return;
        lineY = top - line * ystep;
        UVector result(m_rel->evaluateVector(width));
        if (cancelled()) return;
        typedef QRgb* __restrict QRgbR;
        QRgbR p = reinterpret_cast<QRgb*>(img.scanLine(line));
        for (int x = 0; x < width; ++x) {
            *p++ = result[x] ? fill : 0u;
        }
    }
    emit updated(downsample(img));
}

ImplicitGraph::ImplicitGraph(QObject* parent): Graph2D(parent) {
}

void ImplicitGraph::reset(std::unique_ptr<Equation> rel, const Variable& _x, const Variable& _y) {
    eqn = Sub::create(std::move(rel->a), std::move(rel->b))->simplify();
    x = _x;
    y = _y;
    x.id->type = Variable::Id::Vector;
    y.id->type = Variable::Id::Vector;
    dx = eqn->derivative(x)->simplify();
    dy = eqn->derivative(y)->simplify();
}

void ImplicitGraph::resubstitute() {
    x.id->p = m_px.get();
    y.id->p = m_py.get();
}

QImage ImplicitGraph::draw() {
    QImage img(size().width(), size().height(), QImage::Format_ARGB32_Premultiplied);
    img.fill(qRgba(0, 0, 0, 0));
    QPainter painter(&img);
    // Scale up by a factor of `supersample'
    painter.scale(supersample, supersample);
    for (uz i = 0; i < numPts; ++i) {
        painter.fillRect(QRectF(QPointF(m_px[i], m_py[i]) * transform(), QSizeF(0, 0)).adjusted(-0.5, -0.5, 0.5, 0.5), color());
    }
    painter.end();
    return downsample(img);
}

inline qreal lengthSquared(QPointF p) {
    return p.x() * p.x() + p.y() * p.y();
}

void ImplicitGraph::compute() {
    const int width = size().width(), height = size().height();
    m_px.reset(VECTOR_ALLOC(width * height));
    m_py.reset(VECTOR_ALLOC(width * height));
    resubstitute();
    
    if (cancelled()) return;
    
    // Seed the graph with its initial set of points.
    QTransform ti = transform().inverted();
    numPts = 0;
    uz size = 0;
    // For every 3x3 block, pick three points: (0,0), (1,2), (2,1).
    for (int by = 0; by < height; by += 3) {
        for (int bx = 0; bx < width; bx += 3) {
            for (int delta = 0; delta < 3; ++delta) {
                const int dx = ((int[]){0, 1, 2})[delta], dy = ((int[]){0, 2, 1})[delta];
                QPointF pt = QPointF(bx + dx, by + dy) * ti;
                m_px[size] = pt.x();
                m_py[size] = pt.y();
                ++size;
            }
        }
        if (cancelled()) return;
    }
    // Compute initial correction vectors.
    UVector gx(dx->evaluateVector(size));
    if (cancelled()) return;
    UVector gy(dy->evaluateVector(size));
    if (cancelled()) return;
    UVector p (eqn->evaluateVector(size));
    if (cancelled()) return;
    for (uz i = 0; i < size; ++i) {
        QPointF pt(m_px[i], m_py[i]);
        QPointF correction = QPointF(gx[i], gy[i]) * (p[i] / (gx[i] * gx[i] + gy[i] * gy[i]));
        QPointF visualCorrection = (pt * transform()) - ((pt - correction) * transform()); // in pixels
        if (lengthSquared(visualCorrection) < 200) {
            // Take points that are sufficiently close.
            // While we're at it, perform the correction.
            pt -= correction;
            m_px[numPts] = pt.x();
            m_py[numPts] = pt.y();
            ++numPts;
        }
    }
    if (cancelled()) return;
    
    // Produce an image from this first step.
    emit updated(draw());
    
    // Now do it again.
    forever {
        if (cancelled()) return;
        uz size = numPts;
        UVector p_gx (dx->evaluateVector(size));
        if (cancelled()) return;
        UVector p_gy (dy->evaluateVector(size));
        if (cancelled()) return;
        UVector p_p  (eqn->evaluateVector(size));
        if (cancelled()) return;
        const VectorR
            gx = p_gx.get(),
            gy = p_gy.get(),
            p  = p_p.get(),
            px = m_px.get(),
            py = m_py.get();
        VECTOR_LOOP {
            auto _gx = V(gx), _gy = V(gy);
            auto _p = V(p) / (_gx * _gx + _gy * _gy);
            // delta-P = (gx, gy) * p / (gx^2 + gy^2) = (_gx * _p, _gy * _p)
            // update points
            V(px) -= _gx * _p;
            V(py) -= _gy * _p;
        }
        if (cancelled()) return;
        emit updated(draw());
    }
}

BOOST_CONSTEXPR_OR_CONST uz ParametricGraph::numPts;

ParametricGraph::ParametricGraph(QObject* parent): Graph2D(parent), pts(VECTOR_ALLOC(numPts)) {
}

void ParametricGraph::reset(EPtr _x, EPtr _y, const Variable& _t, Number _tMin, Number _tMax) {
    x = _x->simplify();
    y = _y->simplify();
    t = _t;
    t.id->type = Variable::Id::Vector;
    t.id->p = pts.get();
    tMin = _tMin;
    tMax = _tMax;
}

void ParametricGraph::compute() {
    // tMin and tMax might not be given.
    // This occurs when either x or y is equal to t, so we can figure out the bounds by looking at the transform matrix.
    if (tMin == 0 && tMax == 0) {
        Variable* vx = dynamic_cast<Variable*>(x.get());
        Variable* vy = dynamic_cast<Variable*>(y.get());
        QTransform ti = transform().inverted();
        if (vx && *vx == t) {
            // x = t
            tMin = (QPointF(0, 0) * ti).x();
            tMax = (QPointF(size().width()/supersample, 0) * ti).x();
        } else if (vy && *vy == t) {
            // y = t
            tMax = (QPointF(0, 0) * ti).y();
            tMin = (QPointF(0, size().height()/supersample) * ti).y();
        } // Otherwise, this is a bug, and we'll hit the assertion below.
        Q_ASSERT(tMax > tMin);
    }
    if (cancelled()) return;
    
    // Set up.
    std::uniform_real_distribution<Number> distribution(tMin, tMax);
    std::mt19937 engine;
    QImage img(size().width(), size().height(), QImage::Format_ARGB32_Premultiplied);
    img.fill(qRgba(0, 0, 0, 0));
    QPainter painter(&img);
    painter.scale(supersample, supersample);

    // Start iterating.
    forever {
        if (cancelled()) return;
        // Generate some t-values.
        for (uz i = 0; i < numPts; ++i) pts[i] = distribution(engine);
        if (cancelled()) return;
        // Try to evaluate x and y in parallel.
        QFuture<Vector> fy = QtConcurrent::run(y.get(), &Expression::evaluateVector, numPts);
        UVector vx(x->evaluateVector(numPts));
        fy.waitForFinished();
        UVector vy(fy.result());
        if (cancelled()) return;
        // Draw points.
        for (uz i = 0; i < numPts; ++i) {
            painter.fillRect(QRectF(QPointF(vx[i], vy[i]) * transform(), QSizeF(0, 0)).adjusted(-0.5, -0.5, 0.5, 0.5), color());
        }
        if (cancelled()) return;
        emit updated(downsample(img));
    }
}

