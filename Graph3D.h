#ifndef _GRAPH3D_H_
#define _GRAPH3D_H_

#include "Graph.h"

#include "Render3D.h"
#include "align.h"

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

#endif
