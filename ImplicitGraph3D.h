#ifndef _IMPLICITGRAPH3D_H_
#define _IMPLICITGRAPH3D_H_

#include "Graph3D.h"

#include "Expression.h"
#include <QPixmap>
#include <array>

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

#endif