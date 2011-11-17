#ifndef _GRAPHPROPERTIES_H_
#define _GRAPHPROPERTIES_H_

#include "ui_GraphProperties.h"
#include <QGroupBox>
#include "Expression.h"

class GraphProperties : public QGroupBox, private Ui_GraphProperties {
    Q_OBJECT
public:
    enum Type {
        Relation,
        Parametric
    };
private:
    Variable x, y, t;
    Type lastType;
public:
    GraphProperties(QWidget* parent = 0);
signals:
    void equationChanged(QObject* id, Equation* str, Variable x, Variable y);
    void inequalityChanged(QObject* id, Inequality* str, Variable x, Variable y);
    void parametricChanged(QObject* id, Expression* x, Expression* y, Variable t, Number tMin, Number tMax);
//    void deleteRequest();
public slots:
    void textChanged();
    void deletePressed();
    void setErrorMsg(QString str);
};

#endif
