#ifndef _GRAPHPROPERTIES_H_
#define _GRAPHPROPERTIES_H_

#include "ui_GraphProperties.h"
#include <QGroupBox>

class Graph;

class GraphProperties : public QGroupBox, private Ui_GraphProperties {
    Q_OBJECT
public:
    enum Type {
        Relation,
        Parametric
    };
private:
    Type lastType;
    int color;
public:
    GraphProperties(QWidget* parent = 0);
    virtual ~GraphProperties();
signals:
    void graphChanged(QObject* id, Graph* graph);
//    void deleteRequest();
public slots:
    void textChanged();
    void deletePressed();
    void setErrorMsg(QString str);
};

#endif
