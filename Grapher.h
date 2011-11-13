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

class Graph : QObject {
    Q_OBJECT
    Grapher* parent;
    QFutureWatcher<QImage>* watcher;
    int width, height;
public:
    bool valid;
    boost::scoped_ptr<Expression> eqn, _dx, _dy;
    Variable x, y;
    boost::scoped_array<Number> px, py;
    std::size_t numPts;
    boost::scoped_ptr<Expression> sub, dx, dy;
    QImage img;
    QTransform transform;
    QFuture<QImage> future;
    QFuture<void> setupFuture;
    Graph(Grapher* p);
    ~Graph();
    void reset(const Equation& rel, const Variable& x, const Variable& y);
    void restart(const QTransform& t, const QTransform& ti, int width, int height);
    void resubstitute();
    void iterate();
    QImage draw();
public slots:
    void iterateAgain();
};

class Grapher : public QWidget {
	Q_OBJECT
    Variable x, y;
    QRectF sceneRect;
    QTransform transform;
public:
    QMap<QObject*, Graph*> graphs;
    QMap<QObject*, QFuture<void> > graphFutures;
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

