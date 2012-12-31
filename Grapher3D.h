#ifndef _GRAPHER3D_H_
#define _GRAPHER3D_H_

#include <QWidget>
#include <QMap>

#include "Expression.h"
#include "Render3D.h"
#include "align.h"

class QPaintEvent;
class QResizeEvent;
class QMouseEvent;
class Graph3D;
class QTimer;
class QLabel;

class Grapher3D : public QWidget {
    Q_OBJECT
    struct AData {
        Vector3D<float> boxa, boxb, light;
        Transform3D baseTransform, rotation, comb;
    };
    Align<AData> m_a;

    bool needsRedraw;
    QTimer* redrawTimer;
    bool showAxes;
    QLabel* diagnostic;

    QPoint mouse;
    
    QMap<QObject*, Graph3D*> graphs;
    QMap<Graph3D*, Buffer3D*> images;

protected:
    virtual void mousePressEvent(QMouseEvent*);
    virtual void mouseMoveEvent(QMouseEvent*);

public:
    Grapher3D(QWidget* parent = NULL);
    ~Grapher3D();
    void deleteGraph(QObject* id);
    void addGraph(QObject* id);
    void paintEvent(QPaintEvent* event);
    void resizeEvent(QResizeEvent* event);
signals:
    void lightSourceChanged(Vector3D<float> light);
public slots:
    void setShowAxes(bool showAxes);
    void idDeleted(QObject* id);
    void changeGraph(QObject* id, Graph3D* graph);
    void resized();
    void setBox(Vector3D<float> boxa, Vector3D<float> boxb);
    void setLightSource(Vector3D<float> light);
    void graphUpdated(Buffer3D* img);
    void scheduleUpdate(bool now = false);
    void scheduledUpdate();
};

#endif

