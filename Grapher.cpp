#include "Grapher.h"

#include <QPainter>
#include <QDebug>
#include <QTimer>

#include "Graph.h"

Grapher::Grapher(QWidget* parent) : QWidget(parent), needsRedraw(false), redrawTimer(new QTimer(this)) {
    redrawTimer->setInterval(1000 / 15);
    connect(redrawTimer, SIGNAL(timeout()), this, SLOT(scheduledUpdate()));
}

Grapher::~Grapher() {
    foreach (Graph* graph, graphs) delete graph;
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
    scheduleUpdate(true);
}

void Grapher::idDeleted(QObject* id) {
    deleteGraph(id);
}

void Grapher::changeEquation(QObject* id, Equation* eqn, Variable x, Variable y) {
    Graph* g_graph = graphs[id];
    ImplicitGraph* graph = qobject_cast<ImplicitGraph*>(g_graph);
    if (graph == NULL) {
        if (g_graph) {
            g_graph->cancel();
            delete g_graph;
        }
        graphs[id] = graph = new ImplicitGraph(this);
        connect(graph, SIGNAL(updated()), SLOT(scheduleUpdate()));
    }
    graph->reset(*eqn, x, y);
    delete eqn;
    graph->setupRestart(transform, width(), height());
}

void Grapher::changeInequality(QObject* id, Inequality* eqn, Variable x, Variable y) {
    Graph* g_graph = graphs[id];
    InequalityGraph* graph = qobject_cast<InequalityGraph*>(g_graph);
    if (graph == NULL) {
        if (g_graph) {
            g_graph->cancel();
            delete g_graph;
        }
        graphs[id] = graph = new InequalityGraph(this);
        connect(graph, SIGNAL(updated()), SLOT(scheduleUpdate()));
    }
    graph->reset(std::unique_ptr<Inequality>(eqn), x, y);
    graph->setupRestart(transform, width(), height());
}

void Grapher::changeParametric(QObject* id, Expression* x, Expression* y, Variable t, Number tMin, Number tMax) {
    Graph* g_graph = graphs[id];
    ParametricGraph* graph = qobject_cast<ParametricGraph*>(g_graph);
    if (graph == NULL) {
        if (g_graph) {
            g_graph->cancel();
            delete g_graph;
        }
        graphs[id] = graph = new ParametricGraph(this);
        connect(graph, SIGNAL(updated()), SLOT(scheduleUpdate()));
    }
    graph->reset(EPtr(x), EPtr(y), t, tMin, tMax);
    graph->setupRestart(transform, width(), height());
}

void Grapher::scheduleUpdate(bool now) {
    if (now || !redrawTimer->isActive()) {
        update();
        needsRedraw = false;
        redrawTimer->start();
    } else {
        needsRedraw = true;
    }
}

void Grapher::scheduledUpdate() {
    if (needsRedraw) {
        update();
        needsRedraw = false;
    } else {
        redrawTimer->stop();
    }
}

