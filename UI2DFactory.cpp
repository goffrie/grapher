#include "UI2DFactory.h"

UI2DFactory::UI2DFactory(QObject* parent) : UIFactory(parent) {
}

Grapher2D* UI2DFactory::grapher() {
    if (!m_grapher.isNull()) return m_grapher.data();
    m_grapher = new Grapher2D;
    return m_grapher.data();
}

Graph2DProperties* UI2DFactory::graphProperties() {
    Graph2DProperties* newGraph = new Graph2DProperties;
    grapher()->addGraph(newGraph);
    connect(newGraph, SIGNAL(graphChanged(QObject*, Graph2D*)), grapher(), SLOT(changeGraph(QObject*, Graph2D*)));
    return newGraph;
}

Window2DSettings* UI2DFactory::windowSettings() {
    Window2DSettings* settings = new Window2DSettings;
    connect(settings, SIGNAL(rectChanged(QRectF)), grapher(), SLOT(setWindow(QRectF)));
    connect(settings, SIGNAL(showAxes(bool)), grapher(), SLOT(setShowAxes(bool)));
    connect(settings, SIGNAL(showGrid(bool)), grapher(), SLOT(setShowGrid(bool)));
    settings->onChanged();
    return settings;
}


