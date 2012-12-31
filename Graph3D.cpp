#include "Graph3D.h"

Graph3D::Graph3D(QObject* parent) : Graph(parent) {
}

void Graph3D::setupRestart(const Transform3D& t,
        QSize size,
        Vector3D<float> boxa,
        Vector3D<float> boxb,
        Vector3D<float> light) {
    stop();
    m_size = size;
    m_a->boxa = boxa;
    m_a->boxb = boxb;
    m_a->light = light;
    m_a->transform = t;
    findEyeRay();
    if (width() > 0 && height() > 0) {
        restart();
    }
}

void Graph3D::findEyeRay() {
    Transform3D inv = m_a->transform.inverted();
    Vector3D<float> pt1(0.f,0.f,1.f);
    Vector3D<float> pt2(0.f,0.f,0.f);
    pt1 = inv * pt1;
    pt2 = inv * pt2;
    m_a->eyeray = (pt2 - pt1).normalized<6>();
}
