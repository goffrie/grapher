#include "Grapher.h"

#include <iostream>
#include <QPainter>
#include <QDebug>
#include <QtGlobal>
#include <QtConcurrentMap>
#include <QtConcurrentRun>
#include <QFutureWatcher>
#include "Parser.h"
#include <boost/bind.hpp>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_nan.h>

Grapher::Grapher(QWidget* parent) : QWidget(parent) {
}

Grapher::~Grapher() {
    foreach (Graph* graph, graphs) delete graph;
}

Graph::Graph(Grapher* parent): QObject(parent), watcher(new QFutureWatcher<QImage>(this)) {
    connect(watcher, SIGNAL(finished()), this, SLOT(iterateAgain()));
}

Graph::~Graph() {
    cancel();
}

void Graph::cancel() {
    future.cancel();
    future.waitForFinished();
}

void Graph::setupRestart(const QTransform& t, int _width, int _height) {
    cancel();
    width = _width;
    height = _height;
    transform = t;
    future = QtConcurrent::run(boost::bind(&Graph::restart, this));
    watcher->setFuture(future);
}

ImplicitGraph::ImplicitGraph(Grapher* p) : Graph(p) {
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
    }
    boost::scoped_array<Number> gx(dx->evaluateVector(size)),
                                gy(dy->evaluateVector(size)),
                                p(sub->evaluateVector(size));
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
    return draw();
}

void ImplicitGraph::iterateAgain() {
    future.waitForFinished();
    if (future.isCanceled()) return;
    m_img = future.result();
    future = QtConcurrent::run(boost::bind(&ImplicitGraph::iterate, this));
    watcher->setFuture(future);
    qobject_cast<QWidget*>(parent())->update();
}

QImage ImplicitGraph::draw() {
    QImage _img(width, height, QImage::Format_ARGB32_Premultiplied);
    _img.fill(qRgba(0, 0, 0, 0));
    QPainter painter(&_img);
    painter.setPen(Qt::black);
    painter.setTransform(transform, true);
    for (std::size_t i = 0; i < numPts; ++i) {
        painter.drawPoint(QPointF(m_px[i], m_py[i]));
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
    std::size_t size = numPts;
    VectorR gx  = dx->evaluateVector(size),
            gy  = dy->evaluateVector(size),
            p   = sub->evaluateVector(size),
            gx2 = new Number[size],
            gy2 = new Number[size],
            gxy = new Number[size],
            px = m_px.get(),
            py = m_py.get();
    // correcting factor = (gx, gy) * p / (gx * gx + gy * gy)
    // (gx2, gy2) = (gx^2, gy^2) in place
    for (std::size_t i = 0; i < size; ++i) {
        gx2[i] = gx[i] * gx[i];
    }
    for (std::size_t i = 0; i < size; ++i) {
        gy2[i] = gy[i] * gy[i];
    }
    // gxy = gx^2+gy^2
    for (std::size_t i = 0; i < size; ++i) {
        gxy[i] = gx2[i] + gy2[i];
    }
    // replace "p" with p/(gx^2 + gy^2)
    for (std::size_t i = 0; i < size; ++i) {
        p[i] = p[i] / gxy[i];
    }
    // replace (gx, gy) with (gx, gy) * p / (gx^2 + gy^2)
    for (std::size_t i = 0; i < size; ++i) {
        gx[i] = gx[i] * p[i];
    }
    for (std::size_t i = 0; i < size; ++i) {
        gy[i] = gy[i] * p[i];
    }
    // update points
    for (std::size_t i = 0; i < size; ++i) {
        px[i] = px[i] - gx[i];
    }
    for (std::size_t i = 0; i < size; ++i) {
        py[i] = py[i] - gy[i];
    }
    delete[] gx;
    delete[] gy;
    delete[] p;
    delete[] gx2;
    delete[] gy2;
    delete[] gxy;
    return draw();
}

ParametricGraph::ParametricGraph(Grapher* parent): Graph(parent) {

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
        QFuture<Vector> fx = QtConcurrent::run(boost::bind(&Expression::evaluateVector, sx.get(), num));
        QFuture<Vector> fy = QtConcurrent::run(boost::bind(&Expression::evaluateVector, sy.get(), num));
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
    return _img;
}

void ParametricGraph::draw(Vector vx, Vector vy, size_t n) {
    QPainter painter(&_img);
    painter.setPen(Qt::black);
    painter.setTransform(transform, true);
    for (std::size_t i = 0; i < n; ++i) {
        painter.drawPoint(QPointF(vx[i], vy[i]));
    }
    painter.end();
}

QImage ParametricGraph::iterate() {
    VectorR pt = new Number[numPts - 1];
    VectorR nt = new Number[numPts * 2 - 1];
    VectorR nx = new Number[numPts * 2 - 1];
    VectorR ny = new Number[numPts * 2 - 1];
    std::size_t numT = 0, numNT = 0;
    Number xscale = transform.m11(), yscale = transform.m22();
    bool kept = false;
    for (std::size_t i = 1; i < numPts; ++i) {
        Q_ASSERT(!std::isnan(m_pt[i-1]));
        if (std::isnan(m_pt[i])) {
            if (kept) {
                // fill in gaps with NaN to reduce wasted calculations
                nt[numNT] = GSL_NAN;
                Q_ASSERT(gsl_isnan((nx[numNT] = ny[numNT] = GSL_NAN)));
                ++numNT;
            }
            kept = false;
            ++i;
            continue;
        }
        Number dx = (m_vx[i] - m_vx[i-1]) * xscale, dy = (m_vy[i] - m_vy[i-1]) * yscale;
        if (dx*dx + dy*dy > 0.25) {
            // deepen
            if (!kept) {
                nt[numNT] = m_pt[i-1];
                nx[numNT] = m_vx[i-1];
                ny[numNT] = m_vy[i-1];
                ++numNT;
            }
            nt[numNT] = pt[numT] = (m_pt[i] + m_pt[i-1]) / 2;
            Q_ASSERT(gsl_isnan((nx[numNT] = ny[numNT] = GSL_NAN)));
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
                Q_ASSERT(gsl_isnan((nx[numNT] = ny[numNT] = GSL_NAN)));
                ++numNT;
            }
            kept = false;
        }
    }
    qDebug() << "deepened by" << numT << "old" << numPts << "new" << numNT;
    if (numT == 0) {
        // done.
        delete[] pt;
        delete[] nt;
        delete[] nx;
        delete[] ny;
        m_pt.reset();
        m_vx.reset();
        m_vy.reset();
        numPts = 0;
        return _img;
    }
    VectorR vx, vy;
    {
        EPtr sx, sy;
        {
            Expression::Subst s;
            External T(pt);
            s.insert(std::make_pair(t, &T));
            sx = x->substitute(s);
            sy = y->substitute(s);
        }
        QFuture<Vector> fx = QtConcurrent::run(boost::bind(&Expression::evaluateVector, sx.get(), numT));
        QFuture<Vector> fy = QtConcurrent::run(boost::bind(&Expression::evaluateVector, sy.get(), numT));
        fx.waitForFinished();
        fy.waitForFinished();
        vx = fx.result();
        vy = fy.result();
    }
    QFuture<void> drawer = QtConcurrent::run(boost::bind(&ParametricGraph::draw, this, vx, vy, numT));
    for (std::size_t i = 0, j = 0; j < numNT; ++j) {
        if (gsl_isnan(nt[j])) continue;
        if (nt[j] == pt[i]) {
            Q_ASSERT(gsl_isnan(ny[j]));
            Q_ASSERT(gsl_isnan(nx[j]));
            nx[j] = vx[i];
            ny[j] = vy[i];
            ++i;
        }
    }
    numPts = numNT;
    m_vx.reset(nx);
    m_vy.reset(ny);
    m_pt.reset(nt);
    delete[] pt;
    drawer.waitForFinished();
    delete[] vx;
    delete[] vy;
    return _img;
}


void ParametricGraph::iterateAgain() {
    future.waitForFinished();
    if (future.isCanceled()) return;
    m_img = future.result();
    if (numPts != 0) {
        future = QtConcurrent::run(boost::bind(&ParametricGraph::iterate, this));
        watcher->setFuture(future);
    }
    qobject_cast<QWidget*>(parent())->update();
}

void Grapher::addGraph(QObject* id) {
    graphs.insert(id, NULL);
    connect(id, SIGNAL(destroyed(QObject*)), SLOT(idDeleted(QObject*)));
}

void Grapher::setWindow(QRectF rect) {
    sceneRect = rect;
    resized();
}

void Grapher::resized() {
    QTransform t;
    QSizeF scenesize = sceneRect.size();
    t.scale(width() / scenesize.width(), -height() / scenesize.height());
    QPointF tl = sceneRect.topLeft();
    t.translate(-tl.x(), -scenesize.height()-tl.y());
    transform = t;

    foreach (Graph* graph, graphs) {
        if (!graph) continue;
        graph->setupRestart(t, width(), height());
    }
}

void Grapher::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    foreach (Graph* graph, graphs) {
        if (!graph) continue;
        painter.drawImage(0, 0, graph->img());
    }
}

void Grapher::resizeEvent(QResizeEvent*) {
    resized();
}

void Grapher::deleteGraph(QObject* id) {
    QMap<QObject*, Graph*>::iterator it = graphs.find(id);
    Graph* graph = it.value();
    graph->cancel();
    graph->deleteLater();
    graphs.erase(it);
    update();
}

void Grapher::idDeleted(QObject* id) {
    deleteGraph(id);
}

void Grapher::changeEquation(QObject* id, Equation* eqn, Variable x, Variable y) {
    Graph* g_graph = graphs[id];
    ImplicitGraph* graph = qobject_cast<ImplicitGraph*>(g_graph);
    if (graph == NULL) {
        delete g_graph;
        graphs[id] = graph = new ImplicitGraph(this);
    }
    graph->reset(*eqn, x, y);
    delete eqn;
    graph->setupRestart(transform, width(), height());
}

void Grapher::changeParametric(QObject* id, Expression* x, Expression* y, Variable t, Number tMin, Number tMax) {
    Graph* g_graph = graphs[id];
    ParametricGraph* graph = qobject_cast<ParametricGraph*>(g_graph);
    if (graph == NULL) {
        delete g_graph;
        graphs[id] = graph = new ParametricGraph(this);
    }
    graph->reset(EPtr(x), EPtr(y), t, tMin, tMax);
    graph->setupRestart(transform, width(), height());
}

