#include "UI3DFactory.h"

UI3DFactory::UI3DFactory(QObject* parent) : UIFactory(parent) {
}

Grapher3D* UI3DFactory::grapher() {
    if (!m_grapher.isNull()) return m_grapher.data();
    m_grapher = new Grapher3D;
    return m_grapher.data();
}

Graph3DProperties* UI3DFactory::graphProperties() {
    Graph3DProperties* newGraph = new Graph3DProperties;
    grapher()->addGraph(newGraph);
    connect(newGraph, SIGNAL(graphChanged(QObject*, Graph3D*)), grapher(), SLOT(changeGraph(QObject*, Graph3D*)));
    return newGraph;
}

Window3DSettings* UI3DFactory::windowSettings() {
    Window3DSettings* settings = new Window3DSettings;
    connect(settings, SIGNAL(boxChanged(Vector3D<float>,Vector3D<float>)), grapher(), SLOT(setBox(Vector3D<float>,Vector3D<float>)));
    connect(settings, SIGNAL(lightChanged(Vector3D<float>)), grapher(), SLOT(setLightSource(Vector3D<float>)));
    connect(grapher(), SIGNAL(lightSourceChanged(Vector3D<float>)), settings, SLOT(setLight(Vector3D<float>)));
    connect(settings, SIGNAL(showAxes(bool)), grapher(), SLOT(setShowAxes(bool)));
    settings->onChanged();
    return settings;
}


