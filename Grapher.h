#ifndef _GRAPHER_H_
#define _GRAPHER_H_

#include <QWidget>
#include <QMap>
#include <QFuture>
#include "Expression.h"
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

class QPaintEvent;
class QResizeEvent;

class Grapher;

class Graph : public QObject {
    Q_OBJECT
public:
    bool valid;
    Graph(Grapher* parent);
    virtual ~Graph();
    void setupRestart(const QTransform& t, int width, int height);
    QImage img() { return m_img; }
    void cancel();
protected:
    QFuture<QImage> future;
    QFutureWatcher<QImage>* watcher;
    int width, height;
    QTransform transform;
    QImage m_img;
    virtual QImage restart() = 0;
protected slots:
    virtual void iterateAgain() = 0;
};

class ImplicitGraph : public Graph {
    Q_OBJECT
    
    Variable x, y;
    std::unique_ptr<Expression> eqn, _dx, _dy;
    std::unique_ptr<Expression> sub, dx, dy;
    boost::scoped_array<Number> m_px, m_py;
    std::size_t numPts;
    
    QImage iterate();
public:
    ImplicitGraph(Grapher* parent);
    void reset(const Equation& rel, const Variable& x, const Variable& y);
    void resubstitute();
    QImage draw();
protected:
    virtual QImage restart();
protected slots:
    virtual void iterateAgain();
};

class Grapher : public QWidget {
	Q_OBJECT
    Variable x, y;
    QRectF sceneRect;
    QTransform transform;
public:
    QMap<QObject*, Graph*> graphs;
	Grapher(QWidget* parent = NULL);
    ~Grapher();
    Variable getX() const { return x; }
    Variable getY() const { return y; }
    void deleteGraph(QObject* id);
	void addGraph(QObject* id);
    void paintEvent(QPaintEvent* event);
    void resizeEvent(QResizeEvent* event);
public slots:
    void idDeleted(QObject* id);
    void changeEquation(QObject* id, Equation* eqn);
    void resized();
    void setWindow(QRectF window);
};

#endif

