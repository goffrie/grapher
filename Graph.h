#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <QObject>
#include <QColor>
#include <QFuture>
#include <QFutureWatcher>
#include <QImage>
#include <QMutex>
#include <QTransform>
#include <memory>
#include <random>

#include "Expression.h"

class Graph : public QObject {
    Q_OBJECT
public:
    Graph(QObject* parent = 0);
    virtual ~Graph() { }
    void setupRestart(const QTransform& t, int width, int height);
    
    virtual QImage img() = 0;
    virtual void cancel() = 0;
    virtual void setColor(QColor c) { color = c; }
protected:
    virtual void startThread() = 0;
    int width, height;
    QTransform transform;
    QColor color;
signals:
    void updated();
};

class InequalityGraph : public Graph {
    Q_OBJECT
public:
    InequalityGraph(QObject* parent = 0);
    void reset(std::unique_ptr<Inequality> rel, const Variable& x, const Variable& y);
    virtual QImage img();
    virtual void cancel();
protected:
    virtual void startThread();

    std::unique_ptr<Inequality> rel;
    Variable x, y;
    
    QFuture<void> future;
    QImage m_img;
    QMutex img_mutex;
    bool cancelled;
    void restart();
};

class IteratingGraph : public Graph {
    Q_OBJECT
public:
    IteratingGraph(QObject* parent = 0);
    QImage img() { return m_img; }
    void cancel();
protected:
    virtual void startThread();

    QFuture<QImage> future;
    QFutureWatcher<QImage>* watcher;
    QImage m_img;
    bool cancelled;
    virtual QImage restart() = 0;
protected slots:
    virtual void iterateAgain() = 0;
};

class ImplicitGraph : public IteratingGraph {
    Q_OBJECT
    
    Variable x, y;
    std::unique_ptr<Expression> eqn, _dx, _dy;
    std::unique_ptr<Expression> sub, dx, dy;
    UVector m_px, m_py;
    std::size_t numPts;
    
    QImage iterate();
    void resubstitute();
    QImage draw();
public:
    ImplicitGraph(QObject* parent = 0);
    void reset(std::unique_ptr<Equation> rel, const Variable& x, const Variable& y);
protected:
    virtual QImage restart();
protected slots:
    virtual void iterateAgain();
};

class ParametricGraph : public IteratingGraph {
    Q_OBJECT
    
    Variable t;
    Number tMin, tMax;
    std::unique_ptr<Expression> x, y, sx, sy;
    std::uniform_real_distribution<Number> distribution;
    std::mt19937 engine;
    UVector pts;
    constexpr static std::size_t numPts = 16384;
    QImage _img;

    void draw(Vector vx, Vector vy);
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
