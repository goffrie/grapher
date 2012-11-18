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

class Graph3D : public Graph {
    Q_OBJECT
protected:
    struct AData {
        Transform3D transform;
        Vector3D boxa, boxb, light;
        Vector3D eyeray;
    };
    Buffer3D m_buf;
    Align<AData> m_a;
public:
    Graph3D(QObject* parent = 0);
    void setupRestart(const Transform3D& t,
            int width, int height,
            Vector3D boxa, Vector3D boxb,
            Vector3D light);
    void findEyeRay();
    virtual const Buffer3D& buf() { return m_buf; }
};

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
    WVectorFunc rayfunc_v, d_rayfunc_v;

    UVector line;

    QFuture<void> future;
    bool cancelled;
public:
    ImplicitGraph3D(QObject* parent = 0);
    ~ImplicitGraph3D();
    void reset(std::unique_ptr<Equation> rel,
            const Variable& x, const Variable& y, const Variable& z);
    QPixmap diagnostics(const Transform3D& inv, int px, int py, QSize size);
protected:
    virtual void cancel();
    virtual void startThread();
    void restart();

    void renderLine(const Transform3D& inv, int y);
    bool renderPoint(const Transform3D& inv,
            int px, int py,
            Vector ox, Vector oy, Vector oz);
};

// A parametric graph in three dimensions with two variables.
// Repeatedly plots points by sampling the parameter space.
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

    void restart();
    QFuture<void> future;
    bool cancelled;
public:
    ParametricGraph3D(QObject* parent = 0);
    void reset(std::unique_ptr<Expression> x,
            std::unique_ptr<Expression> y,
            std::unique_ptr<Expression> z,
            const Variable& t, const Variable& u,
            Number tMin, Number tMax,
            Number uMin, Number uMax);
protected:
    virtual void cancel();
    virtual void startThread();
};

#endif
