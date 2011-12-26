#include "Grapher.h"

#include <QPainter>
#include <QDebug>
#include <QTimer>

#include "Graph.h"
#include "util.h"

#include <mmintrin.h>
#include <emmintrin.h>
#include <xmmintrin.h>
typedef __m128i v4si;

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

inline qreal roundOneDigit(qreal n) {
    if (n < 0) return -roundOneDigit(-n);
    if (n == 0) return 0;
    const qreal l = std::log10(n);
    const qreal f = std::pow(10, l - std::floor(l));
    return rnd_d(f) * std::pow(10, std::floor(l));
}

/*QString textRoundOneDigit(qreal n) {
    if (n < 0) return QString("-") + textRoundOneDigit(-n);
    if (n == 0) return QString("0");
    const qreal l = std::log10(n);
    const qreal fl = std::floor(l);
    long il = (long)fl;
    int digit = qRound(std::pow(10, l - fl));
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
}*/

static uint32_t table[65536];
static bool tableGenerated = false;
void genTable() {
    if (!tableGenerated) {
        for (uint32_t b = 1; b < 65536; ++b) {
            table[b] = (((uint64_t) 1 << (uint64_t) 32) + b - 1) / b;
        }
        tableGenerated = true;
    }
}

const v4si zeroq = {0,0};

QImage combine(QList<QImage> images) {
    if (images.length() == 1) return images.first();
    genTable();
    QSize size = images[0].size();
    const int width = size.width();
    const int height = size.height();
    typedef v4si* dataPtr;
    dataPtr data = (dataPtr) aligned_malloc(sizeof(v4si) * width * height);
    memset(data, 0, sizeof(v4si) * width * height);
    foreach (const QImage& img, images) {
        Q_ASSERT(img.format() == QImage::Format_ARGB32_Premultiplied);
        Q_ASSERT(img.size() == size);
        dataPtr __restrict p = data;
        for (int y = 0; y < height; ++y) {
            const uint32_t* __restrict q = reinterpret_cast<const uint32_t*>(img.scanLine(y));
            for (const dataPtr end = p + width; p != end; ++p, ++q) {
                *p = _mm_add_epi32(*p, _mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(*q), zeroq), zeroq));
            }
        }
    }
    uint32_t* __restrict p = reinterpret_cast<uint32_t*>(data);
    QImage ret(size, QImage::Format_ARGB32_Premultiplied);
    for (int y = 0; y < height; ++y) {
        uchar* __restrict q = ret.scanLine(y);
        for (int x = 0; x < width; ++x) {
            const uint32_t multiplier = table[std::max(p[3], uint32_t(255u))];
#define PROCESS_PIXEL *q++ = (uchar) ((uint64_t(*p++ * 255u) * multiplier) >> 32);
            PROCESS_PIXEL PROCESS_PIXEL PROCESS_PIXEL
#undef PROCESS_PIXEL
            *q++ = (uchar) std::min(*p++, uint32_t(255u));
        }
    }
    aligned_free(data);
    return ret;
}

void Grapher::paintEvent(QPaintEvent*) {
    if (width() == 0) return;
    QPainter painter(this);
    painter.fillRect(0, 0, width(), height(), Qt::white);

    QList<QImage> images;
    foreach (Graph* graph, graphs) {
        if (!graph) continue;
        QImage img = graph->img();
        if (img.isNull()) continue;
        if (img.size() != size()) {
            img = img.scaled(size());
        }
        images.append(img);
    }
    if (!images.empty()) painter.drawImage(0, 0, combine(images));

    if (showAxes) {
        painter.setPen(QColor(0, 0, 0, 192));
        const qreal xScale = qreal(sceneRect.right() - sceneRect.left()) / width();
        const qreal yScale = qreal(sceneRect.bottom() - sceneRect.top()) / height();
        const qreal gridPixels = 80;
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
            painter.setPen(QColor(0, 0, 0, 96));
            painter.drawPath(grid);
        }
        painter.setPen(QColor(0, 0, 0, 192));
        painter.drawPath(axes);
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

