#ifndef _UI2DFACTORY_H_
#define _UI2DFACTORY_H_

#include "UIFactory.h"

#include "Grapher2D.h"
#include "Graph2DProperties.h"
#include "Window2DSettings.h"

#include <QPointer>

class UI2DFactory : public UIFactory {
    Q_OBJECT
    QPointer<Grapher2D> m_grapher;
public:
    explicit UI2DFactory(QObject* parent = 0);
    virtual ~UI2DFactory() { }
    virtual Window2DSettings* windowSettings();
    virtual Graph2DProperties* graphProperties();
    virtual Grapher2D* grapher();
};

#endif
