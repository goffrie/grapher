#include "Grapher.h"

#include <QPainter>
#include <QDebug>
#include <QTimer>

#include "Graph.h"
#include "util.h"

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
    update();
}

void Grapher::resized() {
    QTransform t;
    QSizeF scenesize = sceneRect.size();
    t.scale(width() / scenesize.width(), -height() / scenesize.height());
    QPointF tl = sceneRect.bottomLeft();
    t.translate(-tl.x(), -tl.y());
    transform = t;

    foreach (Graph* graph, graphs) {
        if (!graph) continue;
        graph->setupRestart(t, width(), height());
    }
}

qreal roundOneDigit(qreal n) {
    if (n < 0) return -roundOneDigit(-n);
    if (n == 0) return 0;
    const qreal l = std::log10(n);
    const qreal f = std::pow(10, l - std::floor(l));
    return rnd_d(f) * std::pow(10, std::floor(l));
}
QString textRoundOneDigit(qreal n) {
    if (n < 0) return QString("-") + textRoundOneDigit(-n);
    if (n == 0) return QString("0");
    const qreal l = std::log10(n);
    const qreal fl = std::floor(l);
    long il = (long)fl;
    int digit = rnd(std::pow(10, l - fl));
    if (digit > 9) {
        digit /= 10;
        ++il;
    }
    const QString f = QString::number(digit);
    if (il >= 0) {
        return f + QString(il, '0');
    } else {
        return QString("0.") + QString(-il-1, '0') + f;
    }
}

void Grapher::paintEvent(QPaintEvent*) {
    if (width() == 0) return;
    QPainter painter(this);
    painter.fillRect(0, 0, width(), height(), Qt::white);
    
    painter.setPen(Qt::darkGray);
    const qreal xScale = qreal(sceneRect.right() - sceneRect.left()) / width();
    const qreal yScale = qreal(sceneRect.bottom() - sceneRect.top()) / height();
    const qreal gridPixels = 50;
    const QPointF origin = QPointF(0, 0) * transform;
    // draw X axis
    {
        qreal period = roundOneDigit(xScale * gridPixels);
        qreal y = origin.y();
        if (y < 0) y = 0;
        if (y >= height()) y = height()-1;
        painter.drawLine(0, y, width(), y);
        qreal start = std::floor(sceneRect.left() / period);
        qreal end = std::ceil(sceneRect.right() / period);
        int number = rnd_d(end - start) + 1;
        for (int i = 0; i < number; ++i) {
            Number _x = (start + i) * period;
            qreal x = (QPointF(_x, 0) * transform).x();
            painter.drawLine(x, y-3, x, y+3);
            QString text = QString::number(_x);
            if (y*2 >= height()) {
                // top
                painter.drawText(QRectF(x, y-5, 0, 0), Qt::AlignBottom | Qt::AlignHCenter | Qt::TextDontClip, text);
            } else {
                painter.drawText(QRectF(x, y+5, 0, 0), Qt::AlignTop | Qt::AlignHCenter | Qt::TextDontClip, text);
            }
        }
    }
    // draw Y axis
    {
        qreal period = roundOneDigit(yScale * gridPixels);
        qreal x = origin.x();
        if (x < 0) x = 0;
        if (x >= width()) x = width()-1;
        painter.drawLine(x, 0, x, height());
        qreal start = std::floor(sceneRect.top() / period);
        qreal end = std::ceil(sceneRect.bottom() / period);
        int number = rnd_d(end - start) + 1;
        for (int i = 0; i < number; ++i) {
            Number _y = (start + i) * period;
            qreal y = (QPointF(0, _y) * transform).y();
            painter.drawLine(x-3, y, x+3, y);
            QString text = QString::number(_y);
            if (x*2 >= width()) {
                // left side
                painter.drawText(QRectF(x-5, y, 0, 0), Qt::AlignRight | Qt::AlignVCenter | Qt::TextDontClip, text);
            } else {
                painter.drawText(QRectF(x+5, y, 0, 0), Qt::AlignLeft | Qt::AlignVCenter | Qt::TextDontClip, text);
            }
        }
    }
    
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
    if (graph) {
        graph->cancel();
        delete graph;
    }
    graphs.erase(it);
    scheduleUpdate(true);
}

void Grapher::idDeleted(QObject* id) {
    deleteGraph(id);
}

void Grapher::changeGraph(QObject* id, Graph* graph) {
    Graph* g_graph = graphs[id];
    if (g_graph) {
        g_graph->cancel();
        delete g_graph;
    }
    graphs[id] = graph;
    graph->setParent(this);
    connect(graph, SIGNAL(updated()), SLOT(scheduleUpdate()));
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

