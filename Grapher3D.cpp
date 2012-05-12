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
    mouse = event->pos();
    if (diagnostic) {
        foreach (Graph3D* graph, graphs) {
            if (!graph) continue;
            ImplicitGraph3D* g = dynamic_cast<ImplicitGraph3D*>(graph);
            if (!g) continue;
            diagnostic->setPixmap(g->diagnostics(m_a->comb.inverted(), event->x(), event->y(), diagnostic->size()));
            diagnostic->show();
            return;
        }
    }
    QWidget::mousePressEvent(event);
}

void Grapher3D::mouseMoveEvent(QMouseEvent* event) {
    QPoint d = event->pos() - mouse;
    mouse = event->pos();
    AData& a = *m_a;
    Transform3D left = (a.rotation * a.baseTransform);
    a.rotation = Transform3D::rotatorY(d.x() * 0.01f) * Transform3D::rotatorX(d.y() * -0.01f) * a.rotation;
    Vector3D mid = (a.boxa + a.boxb) * 0.5f;
    a.light = (a.rotation * a.baseTransform).inverted() * (left * (a.light - mid)) + mid;
    emit lightSourceChanged(a.light);
    resized();
    
    update();
}

void Grapher3D::addGraph(QObject* id) {
    graphs.insert(id, NULL);
    connect(id, SIGNAL(destroyed(QObject*)), SLOT(idDeleted(QObject*)));
}

void Grapher3D::setBox(Vector3D _boxa, Vector3D _boxb) {
    if (m_a->boxa != _boxa || m_a->boxb != _boxb) {
        m_a->boxa = _boxa;
        m_a->boxb = _boxb;
        resized();
        update();
    }
}

void Grapher3D::setLightSource(Vector3D _light) {
    if (m_a->light != _light) {
        m_a->light = _light;
        resized();
        update();
    }
}


void Grapher3D::resized() {
    Vector3D mid = (m_a->boxa + m_a->boxb) * 0.5f;
    m_a->baseTransform = Transform3D::translator(-mid.x(), -mid.y(), -mid.z()) * Transform3D::isometricTransform;
    m_a->comb = (m_a->rotation * m_a->baseTransform).fit(width(), height(), m_a->boxa.x(), m_a->boxb.x(), m_a->boxa.y(), m_a->boxb.y(), m_a->boxa.z(), m_a->boxb.z());

    foreach (Graph3D* graph, graphs) {
        if (!graph) continue;
        graph->setupRestart(m_a->comb, width(), height(), m_a->boxa, m_a->boxb, m_a->light);
    }
}

void Grapher3D::setShowAxes(bool _showAxes) {
    showAxes = _showAxes;
    update();
}

void Grapher3D::paintEvent(QPaintEvent*) {
    if (width() == 0) return;
    QPainter painter(this);
    Buffer3D buf(width(), height(), m_a->comb);
    if (showAxes) for (int i = 0; i < 3; ++i) for (int j = 0; j < 2; ++j) for (int k = 0; k < 2; ++k) {
        Vector3D a;
        for (int l = 0, n = j; l < 3; ++l) {
            if (l == i) continue;
            a.v[l] = (n?m_a->boxa:m_a->boxb).v[l];
            n = k;
        }
        Vector3D b = a;
        a.v[i] = m_a->boxa.v[i];
        b.v[i] = m_a->boxb.v[i];
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
    graph->setupRestart(m_a->comb, width(), height(), m_a->boxa, m_a->boxb, m_a->light);
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

