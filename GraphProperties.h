#ifndef _GRAPHPROPERTIES_H_
#define _GRAPHPROPERTIES_H_

#include "ui_GraphProperties.h"
#include <QGroupBox>
#include "Expression.h"

class GraphProperties : public QGroupBox, private Ui_GraphProperties {
    Q_OBJECT
    Variable x, y;
public:
    GraphProperties(Variable x, Variable y, QWidget* parent = 0);
signals:
    void equationChanged(QObject* id, Equation* str);
//    void deleteRequest();
public slots:
    void textChanged();
    void deletePressed();
    void setErrorMsg(QString str);
};

#endif
