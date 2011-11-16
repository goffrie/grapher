#include "GraphProperties.h"
#include "Parser.h"
#include "dynamic_unique_cast.h"
#include <QDebug>

GraphProperties::GraphProperties(QWidget* parent): QGroupBox(parent), x("x"), y("y"), t("t") {
    setupUi(this);
    par_tMin->setValidator(new QDoubleValidator());
    par_tMax->setValidator(new QDoubleValidator());
}

void GraphProperties::textChanged() {
    try {
        QWidget* current = graphTypeSwitcher->currentWidget();
        if (current == relation_tab) {
            boost::unordered_map<std::string, Expression*> vars;
            vars.insert(std::make_pair<std::string, Expression*>("x", &x));
            vars.insert(std::make_pair<std::string, Expression*>("y", &y));
            auto eqn = dynamic_unique_cast<Equation>(Parser::parse(rel_equation->text().toStdString(), vars));
            if (!eqn) {
                setErrorMsg(QLatin1String("Didn't get an equation!"));
            } else {
                setErrorMsg(QString());
                emit equationChanged(this, eqn.release(), x, y);
            }
        } else if (current == parametric_tab) {
            double tMin = par_tMin->text().toDouble(), tMax = par_tMax->text().toDouble();
            if (tMax <= tMin || !par_tMin->text().size() || !par_tMax->text().size()) {
                setErrorMsg(QLatin1String("Bad interval!"));
            } else {
                boost::unordered_map<std::string, Expression*> vars;
                vars.insert(std::make_pair<std::string, Expression*>("t", &t));
                auto ex = dynamic_unique_cast<Expression>(Parser::parse(par_x->text().toStdString(), vars));
                if (!ex) {
                    setErrorMsg(QLatin1String("x(t) isn't a function!"));
                } else {
                    auto ey = dynamic_unique_cast<Expression>(Parser::parse(par_y->text().toStdString(), vars));
                    if (!ey) {
                        setErrorMsg(QLatin1String("y(t) isn't a function!"));
                    } else {
                        setErrorMsg(QString());
                        emit parametricChanged(this, ex.release(), ey.release(), t, tMin, tMax);
                    }
                }
            }
        } else {
            qDebug() << "Unknown tab?";
        }
    } catch (const std::exception& e) {
        setErrorMsg(QLatin1String(e.what()));
    }
}

void GraphProperties::deletePressed() {
    deleteLater();
}

void GraphProperties::setErrorMsg(QString str) {
    errorMsg->setText(str);
}
