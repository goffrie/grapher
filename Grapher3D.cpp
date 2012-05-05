#include "Grapher3D.h"

#include <QCoreApplication>
#include <QDebug>
#include <QLabel>
#include <QMouseEvent>
#include <QPainter>
#include <QTimer>

#include "Graph3D.h"
#include "Render3D.h"
#include "util.h"

#include <mmintrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>
typedef __m128i v4si;

Grapher3D::Grapher3D(QWidget* parent) : QWidget(parent), needsRedraw(false), redrawTimer(new QTimer(this)), showAxes(true), diagnostic(nullptr) {
    redrawTimer->setInterval(1000 / 15);
    connect(redrawTimer, SIGNAL(timeout()), this, SLOT(scheduledUpdate()));
    if (QCoreApplication::instance()->arguments().contains("--debug")) {
        diagnostic = new QLabel();
        diagnostic->setMinimumSize(QSize(400, 400));
        diagnostic->setMaximumSize(diagnostic->minimumSize());
    }
}

Grapher3D::~Grapher3D() {
}

void Grapher3D::mousePressEvent(QMouseEvent* event) {
    if (diagnostic) {
        foreach (Graph3D* graph, graphs) {
            if (!graph) continue;
            ImplicitGraph3D* g = dynamic_cast<ImplicitGraph3D*>(graph);
            if (!g) continue;
            diagnostic->setPixmap(g->diagnostics(comb.inverted(), event->x(), event->y(), diagnostic->size()));
            diagnostic->show();
            return;
        }
    }
    mouse = event->pos();
    QWidget::mousePressEvent(event);
}

void Grapher3D::mouseMoveEvent(QMouseEvent* event) {
    QPoint d = event->pos() - mouse;
    mouse = event->pos();
    Transform3D left = (rotation * baseTransform);
    rotation = Transform3D::rotatorY(d.x() * 0.01f) * Transform3D::rotatorX(d.y() * -0.01f) * rotation;
    Vector3D mid = (boxa + boxb) * 0.5f;
    light = (rotation * baseTransform).inverted() * (left * (light - mid)) + mid;
    resized();
    
    update();
}

void Grapher3D::addGraph(QObject* id) {
    graphs.insert(id, NULL);
    connect(id, SIGNAL(destroyed(QObject*)), SLOT(idDeleted(QObject*)));
}

void Grapher3D::setBox(Vector3D _boxa, Vector3D _boxb) {
    if (boxa != _boxa || boxb != _boxb) {
        boxa = _boxa;
        boxb = _boxb;
        resized();
        update();
    }
}

void Grapher3D::setLightSource(Vector3D _light) {
    if (light != _light) {
        light = _light;
        resized();
        update();
    }
}


void Grapher3D::resized() {
    Vector3D mid = (boxa + boxb) * 0.5f;
    baseTransform = Transform3D::translator(-mid.x(), -mid.y(), -mid.z()) * Transform3D::isometricTransform;
    comb = (rotation * baseTransform).fit(width(), height(), boxa.x(), boxb.x(), boxa.y(), boxb.y(), boxa.z(), boxb.z());

    foreach (Graph3D* graph, graphs) {
        if (!graph) continue;
        graph->setupRestart(comb, width(), height(), boxa, boxb, light);
    }
}

void Grapher3D::setShowAxes(bool _showAxes) {
    showAxes = _showAxes;
    update();
}

void Grapher3D::paintEvent(QPaintEvent*) {
    if (width() == 0) return;
    QPainter painter(this);
    Buffer3D buf(width(), height(), comb);
    if (showAxes) for (int i = 0; i < 3; ++i) for (int j = 0; j < 2; ++j) for (int k = 0; k < 2; ++k) {
        Vector3D a;
        for (int l = 0, n = j; l < 3; ++l) {
            if (l == i) continue;
            a.v[l] = (n?boxa:boxb).v[l];
            n = k;
        }
        Vector3D b = a;
        a.v[i] = boxa.v[i];
        b.v[i] = boxb.v[i];
        buf.drawTransformLine(a, b, Qt::black);
    }
    foreach (Graph3D* graph, graphs) {
        if (!graph) continue;
        buf.drawBuffer(0, 0, graph->buf());
    }
    painter.drawImage(0, 0, buf.image());
}

void Grapher3D::resizeEvent(QResizeEvent*) {
    resized();
}

void Grapher3D::deleteGraph(QObject* id) {
    QMap<QObject*, Graph3D*>::iterator it = graphs.find(id);
    Graph3D* graph = it.value();
    if (graph) {
        graph->cancel();
        delete graph;
    }
    graphs.erase(it);
    scheduleUpdate(true);
}

void Grapher3D::idDeleted(QObject* id) {
    deleteGraph(id);
}

void Grapher3D::changeGraph(QObject* id, Graph3D* graph) {
    Graph3D* g_graph = graphs[id];
    if (g_graph) {
        g_graph->cancel();
        delete g_graph;
    }
    graphs[id] = graph;
    graph->setParent(this);
    connect(graph, SIGNAL(updated()), SLOT(scheduleUpdate()));
    graph->setupRestart(comb, width(), height(), boxa, boxb, light);
}

void Grapher3D::scheduleUpdate(bool now) {
    if (now || !redrawTimer->isActive()) {
        update();
        needsRedraw = false;
        redrawTimer->start();
    } else {
        needsRedraw = true;
    }
}

void Grapher3D::scheduledUpdate() {
    if (needsRedraw) {
        update();
        needsRedraw = false;
    } else {
        redrawTimer->stop();
    }
}

