#ifndef _GRAPHER2D_H_
#define _GRAPHER2D_H_

#include <QWidget>
#include <QMap>
#include "Expression.h"

class QPaintEvent;
class QResizeEvent;
class Graph2D;
class QTimer;

class Grapher2D : public QWidget {
    Q_OBJECT
    QRectF sceneRect;
    QTransform transform;
    bool needsRedraw;
    QTimer* redrawTimer;
    bool showAxes;
    bool showGrid;
public:
    QMap<QObject*, Graph2D*> graphs;
    Grapher2D(QWidget* parent = NULL);
    ~Grapher2D();
    void deleteGraph(QObject* id);
    void addGraph(QObject* id);
    void paintEvent(QPaintEvent* event);
    void resizeEvent(QResizeEvent* event);
public slots:
    void setShowAxes(bool showAxes);
    void setShowGrid(bool showGrid);
    void idDeleted(QObject* id);
    void changeGraph(QObject* id, Graph2D* graph);
    void resized();
    void setWindow(QRectF window);
    void scheduleUpdate(bool now = false);
    void scheduledUpdate();
};

#endif

