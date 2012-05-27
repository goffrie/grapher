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
    void setupRestart(const Transform3D& t, int width, int height, Vector3D boxa, Vector3D boxb, Vector3D light);
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
    void reset(std::unique_ptr<Equation> rel, const Variable& x, const Variable& y, const Variable& z);
    QPixmap diagnostics(const Transform3D& inv, int px, int py, QSize size);
protected:
    virtual void cancel();
    virtual void startThread();
    void restart();

    void renderLine(const Transform3D& inv, int y);
    bool renderPoint(const Transform3D& inv, int px, int py, Vector ox, Vector oy, Vector oz);
};

#endif
