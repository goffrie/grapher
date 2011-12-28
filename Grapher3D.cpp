#include "Grapher3D.h"

#include <QPainter>
#include <QDebug>
#include <QTimer>

#include "Graph3D.h"
#include "Render3D.h"
#include "util.h"

#include <mmintrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>
typedef __m128i v4si;

Grapher3D::Grapher3D(QWidget* parent) : QWidget(parent), needsRedraw(false), redrawTimer(new QTimer(this)), showAxes(true) {
    redrawTimer->setInterval(1000 / 15);
    connect(redrawTimer, SIGNAL(timeout()), this, SLOT(scheduledUpdate()));
}

Grapher3D::~Grapher3D() {
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
    transform = isometricTransform(width(), height(), boxa.x(), boxb.x(), boxa.y(), boxb.y(), boxa.z(), boxb.z());

    foreach (Graph3D* graph, graphs) {
        if (!graph) continue;
        graph->setupRestart(transform, width(), height(), boxa, boxb, light);
    }
}

void Grapher3D::setShowAxes(bool _showAxes) {
    showAxes = _showAxes;
    update();
}

void Grapher3D::paintEvent(QPaintEvent*) {
    if (width() == 0) return;
    QPainter painter(this);
    foreach (Graph3D* graph, graphs) {
        if (!graph) continue;
        painter.drawImage(0, 0, graph->buf().image());
        break;
    }
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
    graph->setupRestart(transform, width(), height(), boxa, boxb, light);
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

