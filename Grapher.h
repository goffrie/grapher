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
    
    QFuture<QImage> future;
    QFutureWatcher<QImage>* watcher;
    int width, height;
    Variable x, y;
    boost::scoped_ptr<Expression> eqn, _dx, _dy;
    boost::scoped_ptr<Expression> sub, dx, dy;
    boost::scoped_array<Number> m_px, m_py;
    std::size_t numPts;
    QTransform transform;
    
    QImage restart();
    QImage iterate();
private slots:
    void iterateAgain();
public:
    bool valid;
    QImage img;
    Graph(Grapher* p);
    ~Graph();
    void reset(const Equation& rel, const Variable& x, const Variable& y);
    void setupRestart(const QTransform& t, int width, int height);
    void resubstitute();
    QImage draw();
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

