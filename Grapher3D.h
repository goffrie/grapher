#ifndef _GRAPHER3D_H_
#define _GRAPHER3D_H_

#include <QWidget>
#include <QMap>

#include "Expression.h"
#include "Render3D.h"

class QPaintEvent;
class QResizeEvent;
class QMouseEvent;
class Graph3D;
class QTimer;
class QLabel;

class Grapher3D : public QWidget {
    Q_OBJECT
    Vector3D boxa, boxb, light;
    Transform3D baseTransform, rotation, comb;

    bool needsRedraw;
    QTimer* redrawTimer;
    bool showAxes;
    QLabel* diagnostic;

    QPoint mouse;
protected:
    virtual void mousePressEvent(QMouseEvent*);
    virtual void mouseMoveEvent(QMouseEvent*);
public:
    QMap<QObject*, Graph3D*> graphs;
    Grapher3D(QWidget* parent = NULL);
    ~Grapher3D();
    void deleteGraph(QObject* id);
    void addGraph(QObject* id);
    void paintEvent(QPaintEvent* event);
    void resizeEvent(QResizeEvent* event);
public slots:
    void setShowAxes(bool showAxes);
    void idDeleted(QObject* id);
    void changeGraph(QObject* id, Graph3D* graph);
    void resized();
    void setBox(Vector3D boxa, Vector3D boxb);
    void setLightSource(Vector3D light);
    void scheduleUpdate(bool now = false);
    void scheduledUpdate();
};

#endif

