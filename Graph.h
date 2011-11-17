#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <QObject>
#include <QFuture>
#include <QFutureWatcher>
#include <QImage>
#include <QTransform>
#include <memory>

#include "Expression.h"

class Graph : public QObject {
    Q_OBJECT
public:
    Graph(QObject* parent);
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
    bool cancelled;
    virtual QImage restart() = 0;
protected slots:
    virtual void iterateAgain() = 0;
signals:
    void updated();
};

class ImplicitGraph : public Graph {
    Q_OBJECT
    
    Variable x, y;
    std::unique_ptr<Expression> eqn, _dx, _dy;
    std::unique_ptr<Expression> sub, dx, dy;
    std::unique_ptr<Number[]> m_px, m_py;
    std::size_t numPts;
    
    QImage iterate();
    void resubstitute();
    QImage draw();
public:
    ImplicitGraph(QObject* parent = 0);
    void reset(const Equation& rel, const Variable& x, const Variable& y);
protected:
    virtual QImage restart();
protected slots:
    virtual void iterateAgain();
};

class ParametricGraph : public Graph {
    Q_OBJECT
    
    Variable t;
    Number tMin, tMax;
    std::unique_ptr<Expression> x, y;
    std::unique_ptr<Number[]> m_pt, m_vx, m_vy;
    std::size_t numPts;
    QImage _img;
    
    void draw(Vector vx, Vector vy, std::size_t n);
    QImage iterate();
public:
    ParametricGraph(QObject* parent = 0);
    void reset(std::unique_ptr<Expression> x, std::unique_ptr<Expression> y, const Variable& t, Number tMin, Number tMax);
protected:
    virtual QImage restart();
protected slots:
    virtual void iterateAgain();
};

#endif
