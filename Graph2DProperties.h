#ifndef _GRAPH2DPROPERTIES_H_
#define _GRAPH2DPROPERTIES_H_

#include "ui_Graph2DProperties.h"

class Graph2D;

class Graph2DProperties : public QGroupBox, private Ui_Graph2DProperties {
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
    Graph2DProperties(QWidget* parent = 0);
    virtual ~Graph2DProperties();
signals:
    void graphChanged(QObject* id, Graph2D* graph);
public slots:
    void textChanged();
    void deletePressed();
    void setErrorMsg(QString str);
};

#endif
