#include "Grapher.h"

#include <iostream>
#include <QPainter>
#include <QTimer>
#include <QDebug>
#include <QtConcurrentMap>
#include <QtConcurrentRun>
#include <QFutureWatcher>
#include "Parser.h"
#include <boost/bind.hpp>

Grapher::Grapher(QWidget* parent) : QWidget(parent), x("x"), y("y") {
}

Grapher::~Grapher() {
    foreach (Graph* graph, graphs) delete graph;
}

Graph::Graph(Grapher* p) : parent(p), watcher(new QFutureWatcher<QImage>(this)), valid(false) {
    connect(watcher, SIGNAL(finished()), this, SLOT(iterateAgain()));
}

Graph::~Graph() {
    future.cancel();
    future.waitForFinished();
}

void Graph::reset(const Equation& rel, const Variable& _x, const Variable& _y) {
    future.cancel();
    future.waitForFinished();
    valid = true;
    numPts = 0;
    std::cerr << rel.toString() << std::endl;
    eqn.reset(new Sub(rel.a->copy(), rel.b->copy()));
    std::cerr << eqn->toString() << std::endl;
    x = _x;
    y = _y;
    _dx.reset(eqn->derivative(x));
    _dy.reset(eqn->derivative(y));
    std::cerr << _dx->toString() << '|' << _dy->toString() << std::endl;
    eqn.reset(eqn->simplify());
    _dx.reset(_dx->simplify());
    _dy.reset(_dy->simplify());
    std::cerr << eqn->toString() << std::endl;
    std::cerr << _dx->toString() << '|' << _dy->toString() << std::endl;
    resubstitute();
}

void Graph::setupRestart(const QTransform& t, int _width, int _height) {
    future.cancel();
    future.waitForFinished();
    width = _width;
    height = _height;
    transform = t;
    future = QtConcurrent::run(boost::bind(&Graph::restart, this));
    watcher->setFuture(future);
}

QImage Graph::restart() {
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

void Graph::iterateAgain() {
    future.waitForFinished();
    if (future.isCanceled()) return;
    img = future.result();
    future = QtConcurrent::run(boost::bind(&Graph::iterate, this));
    watcher->setFuture(future);
    parent->update();
}

QImage Graph::draw() {
    QImage _img(width, height, QImage::Format_ARGB32_Premultiplied);
    _img.fill(qRgba(0, 0, 0, 0));
    QPainter painter(&_img);
    painter.setPen(Qt::black);
    for (std::size_t i = 0; i < numPts; ++i) {
        painter.drawPoint(QPointF(m_px[i], m_py[i]) * transform);
    }
    painter.end();
    return _img;
}

void Grapher::addGraph(QObject* id) {
    graphs.insert(id, new Graph(this));
    connect(id, SIGNAL(destroyed(QObject*)), SLOT(idDeleted(QObject*)));
}

void Graph::resubstitute() {
    Expression::Subst s;
    External X(m_px.get()), Y(m_py.get());
    s.insert(std::make_pair(x, &X));
    s.insert(std::make_pair(y, &Y));
    dx.reset(_dx->substitute(s));
    dy.reset(_dy->substitute(s));
    sub.reset(eqn->substitute(s));
}

QImage Graph::iterate() {
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
        if (!graph->valid) continue;
        graph->setupRestart(t, width(), height());
    }
}
/*
void Grapher::iterate() {
    foreach (Graph* graph, graphs) {
        if (!graph->valid) continue;
        graph->iterate();
        graph->redraw();
    }
    
    update();
}*/

void Grapher::paintEvent(QPaintEvent*) {
    QPainter painter(this);
    foreach (Graph* graph, graphs) {
        if (!graph->valid) continue;
        painter.drawImage(0, 0, graph->img);
    }
}

void Grapher::resizeEvent(QResizeEvent*) {
    resized();
}

void Grapher::deleteGraph(QObject* id) {
    graphs.erase(graphs.find(id));
    update();
}

void Grapher::idDeleted(QObject* id) {
    deleteGraph(id);
}

void Grapher::changeEquation(QObject* id, Equation* eqn) {
    Graph* graph = graphs[id];
    graph->reset(*eqn, x, y);
    delete eqn;
    graph->setupRestart(transform, width(), height());
}
