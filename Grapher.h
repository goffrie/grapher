#ifndef _GRAPHER_H_
#define _GRAPHER_H_

#include <QWidget>
#include <QMap>
#include "Expression.h"

class QPaintEvent;
class QResizeEvent;
class Graph;
class QTimer;

class Grapher : public QWidget {
	Q_OBJECT
    QRectF sceneRect;
    QTransform transform;
    bool needsRedraw;
    QTimer* redrawTimer;
    bool showAxes;
    bool showGrid;
public:
    QMap<QObject*, Graph*> graphs;
	Grapher(QWidget* parent = NULL);
    ~Grapher();
    void deleteGraph(QObject* id);
	void addGraph(QObject* id);
    void paintEvent(QPaintEvent* event);
    void resizeEvent(QResizeEvent* event);
public slots:
    void setShowAxes(bool showAxes);
    void setShowGrid(bool showGrid);
    void idDeleted(QObject* id);
    void changeGraph(QObject* id, Graph* graph);
    void resized();
    void setWindow(QRectF window);
    void scheduleUpdate(bool now = false);
    void scheduledUpdate();
};

#endif

