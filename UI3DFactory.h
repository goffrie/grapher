#ifndef _UI3DFACTORY_H_
#define _UI3DFACTORY_H_

#include "UIFactory.h"

#include "Grapher3D.h"
#include "Graph3DProperties.h"
#include "Window3DSettings.h"

#include <QPointer>

class UI3DFactory : public UIFactory {
    Q_OBJECT
    QPointer<Grapher3D> m_grapher;
public:
    explicit UI3DFactory(QObject* parent = 0);
    virtual ~UI3DFactory() { }
    virtual Window3DSettings* windowSettings();
    virtual Graph3DProperties* graphProperties();
    virtual Grapher3D* grapher();
};

#endif
