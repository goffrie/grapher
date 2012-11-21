#include "Grapher2D.h"

#include <QPainter>
#include <QDebug>
#include <QTimer>

#include "Graph2D.h"
#include "util.h"
#include "global.h"

#include <Vc/Vc>

Grapher2D::Grapher2D(QWidget* parent) : QWidget(parent), needsRedraw(false), redrawTimer(new QTimer(this)), showAxes(true), showGrid(true) {
    redrawTimer->setInterval(1000 / 15);
    connect(redrawTimer, SIGNAL(timeout()), this, SLOT(scheduledUpdate()));
}

Grapher2D::~Grapher2D() {
    foreach (Graph2D* graph, graphs) delete graph;
}

void Grapher2D::addGraph(QObject* id) {
    graphs.insert(id, NULL);
    connect(id, SIGNAL(destroyed(QObject*)), SLOT(idDeleted(QObject*)));
}

void Grapher2D::setWindow(QRectF rect) {
    sceneRect = rect;
    resized();
    update();
}

void Grapher2D::resized() {
    QTransform t;
    QSizeF scenesize = sceneRect.size();
    t.scale(width() / scenesize.width(), -height() / scenesize.height());
    QPointF tl = sceneRect.bottomLeft();
    t.translate(-tl.x(), -tl.y());
    transform = t;

    foreach (Graph2D* graph, graphs) {
        if (!graph) continue;
        graph->setupRestart(t, width(), height());
    }
}

void Grapher2D::setShowAxes(bool _showAxes) {
    showAxes = _showAxes;
    update();
}

void Grapher2D::setShowGrid(bool _showGrid) {
    showGrid = _showGrid;
    update();
}

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

static const Vc::int_v zeroq(Vc::Zero);

QImage combine(QList<QImage> images) {
    if (images.length() == 1) return images.first();
    genTable();
    QSize size = images[0].size();
    const int width = size.width();
    const int height = size.height();
    QImage ret(size, QImage::Format_ARGB32_Premultiplied);
    for (int y = 0; y < height; ++y) {
        uchar* __restrict q = ret.scanLine(y);
        for (int x = 0; x < width; ++x) {
            u32 sum[4] = {0};
            foreach (const QImage& img, images) {
                QRgb rgb = img.pixel(x, y);
                sum[0] += rgb & 0xFFu;
                sum[1] += (rgb >> 8u) & 0xFFu;
                sum[2] += (rgb >> 16u) & 0xFFu;
                sum[3] += (rgb >> 24u) & 0xFFu;
            }
            const u32 multiplier = table[std::max(sum[3], u32(255u))];
#define PROCESS_PIXEL(i) *q++ = (uchar) ((u64(sum[i] * 255u) * multiplier) >> 32);
            PROCESS_PIXEL(0) PROCESS_PIXEL(1) PROCESS_PIXEL(2)
#undef PROCESS_PIXEL
            *q++ = (u8) std::min(sum[3], u32(255u));
        }
    }
    /*
    typedef Vc::int_v* dataPtr;
    dataPtr data = (dataPtr) aligned_malloc(sizeof(Vc::int_v) * width * height);
    memset(data, 0, sizeof(Vc::int_v) * width * height);
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
    */
    return ret;
}

void Grapher2D::paintEvent(QPaintEvent*) {
    if (width() == 0) return;
    QPainter painter(this);
    painter.fillRect(0, 0, width(), height(), Qt::white);

    QList<QImage> images;
    foreach (Graph2D* graph, graphs) {
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

void Grapher2D::resizeEvent(QResizeEvent*) {
    resized();
}

void Grapher2D::deleteGraph(QObject* id) {
    QMap<QObject*, Graph2D*>::iterator it = graphs.find(id);
    Graph2D* graph = it.value();
    if (graph) {
        graph->cancel();
        delete graph;
    }
    graphs.erase(it);
    scheduleUpdate(true);
}

void Grapher2D::idDeleted(QObject* id) {
    deleteGraph(id);
}

void Grapher2D::changeGraph(QObject* id, Graph2D* graph) {
    Graph2D* g_graph = graphs[id];
    if (g_graph) {
        g_graph->cancel();
        delete g_graph;
    }
    graphs[id] = graph;
    graph->setParent(this);
    connect(graph, SIGNAL(updated()), SLOT(scheduleUpdate()));
    graph->setupRestart(transform, width(), height());
}

void Grapher2D::scheduleUpdate(bool now) {
    if (now || !redrawTimer->isActive()) {
        update();
        needsRedraw = false;
        redrawTimer->start();
    } else {
        needsRedraw = true;
    }
}

void Grapher2D::scheduledUpdate() {
    if (needsRedraw) {
        update();
        needsRedraw = false;
    } else {
        redrawTimer->stop();
    }
}
