#ifndef _GRAPH3D_H_
#define _GRAPH3D_H_

#include "Graph.h"

#include "Expression.h"
#include "Render3D.h"

#include <QColor>
#include <QFuture>

class Graph3D : public Graph {
    Q_OBJECT
protected:
    Transform3D m_transform;
    Buffer3D m_buf;
    Vector3D m_boxa, m_boxb, m_light;
    Vector3D m_eyeray;
public:
    Graph3D(QObject* parent = 0);
    void setupRestart(const Transform3D& t, int width, int height, Vector3D boxa, Vector3D boxb, Vector3D light);
    void findEyeRay();
    virtual const Buffer3D& buf() { return m_buf; }
};

class Sphere : public Graph3D {
    Q_OBJECT
public:
    Sphere(QObject* parent = 0);
    virtual void cancel();
protected:
    virtual void startThread();
};

class ImplicitGraph3D : public Graph3D {
    Q_OBJECT
private:
    Variable x, y, z, tv;
    EPtr func;
    float v1[3], dv[3];
    float te[128] __attribute__((aligned(16)));
    EPtr rayfunc[2];
    EPtr d_rayfunc;
    EPtr dx, dy, dz;

    UVector line;

    QFuture<void> future;
    bool cancelled;
public:
    ImplicitGraph3D(QObject* parent = 0);
    void reset(std::unique_ptr<Equation> rel, const Variable& x, const Variable& y, const Variable& z);
protected:
    virtual void cancel();
    virtual void startThread();
    void restart();

    void renderLine(const Transform3D& inv, int y);
    bool renderPoint(const Transform3D& inv, int py, int px, Vector ox, Vector oy, Vector oz);
    QImage diagnostics(const Transform3D& inv, int py, int px, Vector ox, Vector oy, Vector oz, QSize size);
};

#endif
