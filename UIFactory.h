#ifndef _UIFACTORY_H_
#define _UIFACTORY_H_

#include <QObject>

class QWidget;

class UIFactory : public QObject {
    Q_OBJECT
public:
    explicit UIFactory(QObject* parent = 0) : QObject(parent) { }
    virtual ~UIFactory() { }
    virtual QWidget* windowSettings() = 0;
    virtual QWidget* graphProperties() = 0;
    virtual QWidget* grapher() = 0;
};

#endif
