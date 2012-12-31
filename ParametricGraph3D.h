#ifndef _PARAMETRICGRAPH3D_H_
#define _PARAMETRICGRAPH3D_H_

#include "Graph3D.h"

#include "global.h"

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
    BOOST_STATIC_CONSTEXPR uz numPts = 16384;

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