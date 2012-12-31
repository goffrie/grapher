#ifndef _GRAPH2D_H_
#define _GRAPH2D_H_

#include "Graph.h"
#include "Expression.h"

#include <QFuture>
#include <QFutureWatcher>
#include <QImage>
#include <QMutex>
#include <QTransform>

#include <memory>
#include <random>

#include <boost/config/suffix.hpp>

/**
 * Abstract base class for 2D graphs.
 * Provides an interface for changing the graph settings (e.g. resizing).
 * To subclasses, provides a transformation matrix converting from
 * graph coordinates to image coordinates.
 */
class Graph2D: public Graph {
    Q_OBJECT
    QTransform m_transform;
    QSize m_size;

public:
    Graph2D(QObject* parent = 0);

protected:
    /// Gets the transformation matrix.
    /// Converts from right-handed graph coordinates to left-handed image coordinates.
    QTransform transform() const { return m_transform; }

    /// Gets the width and height of the graph in pixels.
    QSize size() const { return m_size; }

public slots:
    /// Start redrawing the graph with a new view.
    void setupRestart(const QTransform& t, QSize size);

signals:
    /// Indicates that this graph has a new image to draw.
    void updated(QImage);

public:
    /// A knob adjusting how much supersampling should be done.
    /// A value of N means that each pixel averages over an NxN grid of samples.
    BOOST_CONSTEXPR_OR_CONST static int supersample = 2;
};

/**
 * A graph depicting an inequality in 2D Cartesian coordinates.
 */
class InequalityGraph: public Graph2D {
    Q_OBJECT
    std::unique_ptr<Inequality> m_rel;
    Variable x, y;
    QImage m_img;
public:
    InequalityGraph(QObject* parent = 0);

    /// Initializes the graph.
    /// @param[in] rel The inequality to plot.
    /// @param[in] x,y The x and y variables bound in `rel'.
    void reset(std::unique_ptr<Inequality> rel, const Variable& x, const Variable& y);
protected:
    virtual void compute() override;
};

/**
 * A graph of a curve described by an equation in 2D Cartesian coordinates.
 */
class ImplicitGraph: public Graph2D {
    Q_OBJECT

    Variable x, y;
    EPtr eqn, dx, dy;
    UVector m_px, m_py;
    uz numPts;

    /// Rebind \a m_px and \a m_py to \a x and \a y.
    void resubstitute();
    /// Create an image out of the current set of points.
    QImage draw();

public:
    ImplicitGraph(QObject* parent = 0);

    /// Initializes the graph.
    /// @param[in] rel The equation to plot.
    /// @param[in] x,y The x and y variables bound in `rel'.
    void reset(std::unique_ptr<Equation> rel, const Variable& x, const Variable& y);
protected:
    virtual void compute() override;
};

/**
 * A graph of a curve described by parametric equations in a single variable.
 */
class ParametricGraph: public Graph2D {
    Q_OBJECT

    Variable t;
    Number tMin, tMax;
    EPtr x, y;
    UVector pts;
    BOOST_STATIC_CONSTEXPR uz numPts = 16384;

public:
    ParametricGraph(QObject* parent = 0);
    
    /// Initializes the graph.
    /// @param[in] x,y The x and y coordinates of each point, as functions of \a t.
    /// @param[in] t The variable to vary.
    /// @param[in] tMin,tMax The bounds on \a t.
    void reset(EPtr x, EPtr y, const Variable& t, Number tMin, Number tMax);
protected:
    virtual void compute() override;
};

#endif
