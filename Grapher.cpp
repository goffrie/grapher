#include "Grapher.h"

#include <QPainter>
#include <QDebug>
#include <QTimer>

#include "Graph.h"
#include "util.h"

Grapher::Grapher(QWidget* parent) : QWidget(parent), needsRedraw(false), redrawTimer(new QTimer(this)), showAxes(true), showGrid(true) {
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

void Grapher::setShowAxes(bool _showAxes) {
    showAxes = _showAxes;
    update();
}

void Grapher::setShowGrid(bool _showGrid) {
    showGrid = _showGrid;
    update();
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
    
    if (showAxes) {
        painter.setPen(Qt::darkGray);
        const qreal xScale = qreal(sceneRect.right() - sceneRect.left()) / width();
        const qreal yScale = qreal(sceneRect.bottom() - sceneRect.top()) / height();
        const qreal gridPixels = 50;
        const QPointF origin = QPointF(0, 0) * transform;
        QPainterPath grid;
        QPainterPath axes;
        // draw X axis
        {
            qreal period = roundOneDigit(xScale * gridPixels);
            qreal y = origin.y();
            if (y < 0) y = 0;
            if (y >= height()) y = height()-1;
            axes.moveTo(0, y);
            axes.lineTo(width(), y);
            qreal start = std::floor(sceneRect.left() / period);
            qreal end = std::ceil(sceneRect.right() / period);
            int number = rnd_d(end - start) + 1;
            for (int i = 0; i < number; ++i) {
                Number _x = (start + i) * period;
                qreal x = (QPointF(_x, 0) * transform).x();
                if (showGrid) {
                    grid.moveTo(x, 0);
                    grid.lineTo(x, height());
                }
                axes.moveTo(x, y-3);
                axes.lineTo(x, y+3);
                QString text = QString::number(_x);
                int flag = 0;
                if (showGrid)
                    flag |= (x*2 >= width() ? Qt::AlignRight : Qt::AlignLeft);
                else 
                    flag |= Qt::AlignHCenter;
                flag |= Qt::TextDontClip;
                if (y*2 >= height()) {
                    // top
                    painter.drawText(QRectF(x, y-5, 0, 0), Qt::AlignBottom | flag, text);
                } else {
                    painter.drawText(QRectF(x, y+5, 0, 0), Qt::AlignTop | flag, text);
                }
            }
        }
        // draw Y axis
        {
            qreal period = roundOneDigit(yScale * gridPixels);
            qreal x = origin.x();
            if (x < 0) x = 0;
            if (x >= width()) x = width()-1;
            axes.moveTo(x, 0);
            axes.lineTo(x, height());
            qreal start = std::floor(sceneRect.top() / period);
            qreal end = std::ceil(sceneRect.bottom() / period);
            int number = rnd_d(end - start) + 1;
            for (int i = 0; i < number; ++i) {
                Number _y = (start + i) * period;
                qreal y = (QPointF(0, _y) * transform).y();
                if (showGrid) {
                    grid.moveTo(0, y);
                    grid.lineTo(width(), y);
                }
                axes.moveTo(x-3, y);
                axes.lineTo(x+3, y);
                QString text = QString::number(_y);
                int flag = 0;
                if (showGrid)
                    flag |= (y*2 >= height() ? Qt::AlignBottom : Qt::AlignTop);
                else 
                    flag |= Qt::AlignVCenter;
                flag |= Qt::TextDontClip;
                if (x*2 >= width()) {
                    // left side
                    painter.drawText(QRectF(x-5, y, 0, 0), Qt::AlignRight | flag, text);
                } else {
                    painter.drawText(QRectF(x+5, y, 0, 0), Qt::AlignLeft | flag, text);
                }
            }
        }
        if (showGrid) {
            painter.setPen(Qt::lightGray);
            painter.drawPath(grid);
        }
        painter.setPen(Qt::darkGray);
        painter.drawPath(axes);
    }
    
    foreach (Graph* graph, graphs) {
        if (!graph) continue;
        painter.drawImage(QRectF(QPointF(0, 0), size()), graph->img());
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

