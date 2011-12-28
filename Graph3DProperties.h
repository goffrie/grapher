#ifndef _GRAPH3DPROPERTIES_H_
#define _GRAPH3DPROPERTIES_H_

#include "ui_Graph3DProperties.h"

class Graph3D;

class Graph3DProperties : public QGroupBox, private Ui_Graph3DProperties {
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
    Graph3DProperties(QWidget* parent = 0);
    virtual ~Graph3DProperties();
signals:
    void graphChanged(QObject* id, Graph3D* graph);
public slots:
    void textChanged();
    void deletePressed();
    void setErrorMsg(QString str);
};

#endif
