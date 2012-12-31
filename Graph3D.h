#ifndef _GRAPH3D_H_
#define _GRAPH3D_H_

#include "Graph.h"

#include "Expression.h"
#include "Render3D.h"
#include "align.h"

#include <QColor>
#include <QFuture>
#include <QPixmap>

#include <array>

/**
 * Abstract base class for 3D graphs.
 * Provides an interface for changing the graph settings
 * (transform, size, window, lighting).
 * Exposes this information to subclasses.
 */
class Graph3D : public Graph {
    Q_OBJECT
protected:
    struct AData {
        /// Transformation from 3D coordinates to screen coordinates (with Z-component).
        /// Will be an axonometric transformation.
        Transform3D transform;
        /// Window settings.
        Vector3D<float> boxa, boxb;
        /// The normalized direction of the light.
        Vector3D<float> light;
        /// The normalized direction in which the camera is pointing.
        Vector3D<float> eyeray;
    };

    /// Figures out which way the camera is pointing.
    /// Uses the transformation matrix.
    void findEyeRay();

private:
    /// Holds vector quantities, which should be aligned.
    Align<AData> m_a;
    
    /// The screen size of the graph, in pixels.
    QSize m_size;

public:
    Graph3D(QObject* parent = 0);
    
    // Accessors.
    Transform3D transform() const { return m_a->transform; }
    Vector3D<float> boxa() const { return m_a->boxa; }
    Vector3D<float> boxb() const { return m_a->boxb; }
    Vector3D<float> light() const { return m_a->light; }
    Vector3D<float> eyeray() const { return m_a->eyeray; }
    QSize size() const { return m_size; }
    int width() const { return size().width(); }
    int height() const { return size().height(); }

public slots:
    /// Resets the graph with new settings.
    /// @param[in] transform The transformation matrix from scene coordinates to screen coordinates.
    /// @param[in] size The size of the graph in pixels.
    /// @param[in] boxa,boxb The lower and upper x,y,z coordinates of the corners of the view box.
    /// @param[in] light The direction of the lighting (does not need to be normalized).
    void setupRestart(const Transform3D& transform,
            QSize size,
            Vector3D<float> boxa, Vector3D<float> boxb,
            Vector3D<float> light);

signals:
    /// Indicates that this graph has a new image to draw.
    /// Buffer allocated with new.
    void updated(Buffer3D*);
};

/**
 * A graph of a surface described by an equation in 3D Cartesian coordinates.
 */
class ImplicitGraph3D : public Graph3D {
    Q_OBJECT
private:
    Variable x, y, z, tv, v1x, v1y, v1z, dvx, dvy, dvz;
    EPtr func;
    float v1[3], dv[3];
    Align<std::array<float, 256>> te;
    EPtr rayfunc[3];
    std::unique_ptr<Polynomial> polyrayfunc;
    EPtr d_rayfunc;
    EPtr dx, dy, dz;
    std::vector<WEvalFunc> polyrayfunc_e;
    WEvalFunc rayfunc_e[3];

    UVector line;
public:
    ImplicitGraph3D(QObject* parent = 0);
    ~ImplicitGraph3D();

    /// Initializes the graph.
    /// @param[in] rel The 3D equation to plot, in Cartesian coordinates.
    /// @param[in] x,y,z The x, y, and z variables bound in \a rel.
    void reset(std::unique_ptr<Equation> rel,
            const Variable& x, const Variable& y, const Variable& z);
    
    /// Render a diagnostic image for one pixel.
    /// @param[in] size The size of the image to return.
    QPixmap diagnostics(const Transform3D& inv, int px, int py, QSize size);

protected:
    virtual void compute() override;

private:
    /// Render one horizontal line of pixels.
    /// @param[in] inv The matrix inverse of \a transform().
    /// @param[in] y The line to render.
    void renderLine(const Transform3D& inv, int y);
    
    /// Raytrace one pixel.
    /// @param[in] inv The matrix inverse of \a transform().
    /// @param[in] px,py The 2D (screen) coordinates of the pixel.
    /// @param[out] ox,oy,oz The 3D (scene) coordinates of the raytraced point.
    /// @return Returns true if a point was found on this ray, false otherwise.
    bool renderPoint(const Transform3D& inv,
            int px, int py,
            Vector ox, Vector oy, Vector oz);
};

/**
 * A graph of a surface in three dimensions described by parametric equations in two variables.
 * Repeatedly plots points by randomly sampling the parameter space.
 */
class ParametricGraph3D : public Graph3D {
    Q_OBJECT

    // Parameters
    Variable t, u;
    Number tMin, tMax, uMin, uMax;
    // Functions x(t,u), y(t,u), z(t,u)
    EPtr x, y, z;
    // Partial derivatives with respect to t, u
    EPtr x_t, y_t, z_t, x_u, y_u, z_u;
    // RNG
    std::uniform_real_distribution<Number> tDist, uDist;
    std::mt19937 engine;

    // Scratch space to hold `t' and `u' values
    UVector tPts, uPts;
    BOOST_STATIC_CONSTEXPR std::size_t numPts = 16384;

public:
    ParametricGraph3D(QObject* parent = 0);

    /// Initializes the graph.
    /// @param[in] x,y,z The x, y, z coordinates of each point, as functions of \a t and \a u.
    /// @param[in] t,u The parameters.
    /// @param[in] tMin,tMax,uMin,uMax The bounds on \a t and \a u.
    void reset(EPtr x, EPtr y, EPtr z,
            const Variable& t, const Variable& u,
            Number tMin, Number tMax,
            Number uMin, Number uMax);

protected:
    virtual void compute() override;
};

#endif
