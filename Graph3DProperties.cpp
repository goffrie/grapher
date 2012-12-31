#include "Graph3DProperties.h"

#include "Parser.h"
#include "Expression.h"
#include "ImplicitGraph3D.h"
#include "ParametricGraph3D.h"

#include <exception>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_math.h>

#include "dynamic_unique_cast.h"
#include "util.h"

static int nextColor = 0;
static std::set<int> recoveredColors;

static QColor getColor(int c) {
    const static qreal factor = -(3 - std::sqrt(5)) / 2;
    const qreal fc = factor * c;
    return QColor::fromHsvF(fc - std::floor(fc), 0.9, 0.8);
}

Graph3DProperties::Graph3DProperties(QWidget* parent): QGroupBox(parent) {
    setupUi(this);
    par_tMin->setValidator(new QDoubleValidator());
    par_tMax->setValidator(new QDoubleValidator());
    par_uMin->setValidator(new QDoubleValidator());
    par_uMax->setValidator(new QDoubleValidator());
    if (recoveredColors.size() > 0) {
        color = *recoveredColors.begin();
        recoveredColors.erase(recoveredColors.begin());
    } else {
        color = nextColor++;
    }
    colorButton->setIcon(circle(16, getColor(color)));
    QPalette p = colorButton->palette();
    p.setColor(QPalette::Button, getColor(color));
    colorButton->setPalette(p);
}

Graph3DProperties::~Graph3DProperties() {
    recoveredColors.insert(color);
}

struct InvalidInputException : public std::exception {
    std::string str;
    InvalidInputException(std::string s) : str(s) { }
    virtual ~InvalidInputException() throw() { }
    const char* what() const throw() {
        return str.c_str();
    }
};

void Graph3DProperties::textChanged() {
    try {
        std::unordered_map<std::string, Expression*> vars;
        Constant _pi(M_PI);
        Constant _e(M_E);
        vars.insert(std::make_pair("pi", &_pi));
        vars.insert(std::make_pair("e", &_e));
        QWidget* current = graphTypeSwitcher->currentWidget();
        if (current == relation_tab) {
            Variable x("x"), y("y"), z("z");
            vars.insert(std::make_pair<std::string, Expression*>("x", &x));
            vars.insert(std::make_pair<std::string, Expression*>("y", &y));
            vars.insert(std::make_pair<std::string, Expression*>("z", &z));
            std::unique_ptr<Thing> thing(Parser::parse(rel_equation->text().toStdString(), vars));
            std::unique_ptr<Equation> eqn = dynamic_maybe_unique_cast<Equation>(thing);
            if (eqn) {
                setErrorMsg(QString());
                ImplicitGraph3D* graph = new ImplicitGraph3D;
                graph->reset(std::move(eqn), x, y, z);
                graph->setColor(getColor(color));
                emit graphChanged(this, graph);
            } else {
                throw InvalidInputException("Didn't get an equation!");
            }
        } else if (current == parametric_tab) {
            double tMin = par_tMin->text().toDouble(), tMax = par_tMax->text().toDouble();
            if (tMax <= tMin || !par_tMin->text().size() || !par_tMax->text().size()) {
                throw InvalidInputException("Bad interval for t!");
            }
            double uMin = par_uMin->text().toDouble(), uMax = par_uMax->text().toDouble();
            if (uMax <= uMin || !par_uMin->text().size() || !par_uMax->text().size()) {
                throw InvalidInputException("Bad interval for u!");
            }
            Variable t("t"), u("u");
            vars.insert(std::make_pair<std::string, Expression*>("t", &t));
            vars.insert(std::make_pair<std::string, Expression*>("u", &u));
            auto ex = dynamic_unique_cast<Expression>(Parser::parse(par_x->text().toStdString(), vars));
            if (!ex) {
                throw InvalidInputException("x(t,u) isn't a function!");
            }
            auto ey = dynamic_unique_cast<Expression>(Parser::parse(par_y->text().toStdString(), vars));
            if (!ey) {
                throw InvalidInputException("y(t,u) isn't a function!");
            }
            auto ez = dynamic_unique_cast<Expression>(Parser::parse(par_z->text().toStdString(), vars));
            if (!ez) {
                throw InvalidInputException("z(t,u) isn't a function!");
            }
            setErrorMsg(QString());
            ParametricGraph3D* graph = new ParametricGraph3D;
            graph->reset(std::move(ex), std::move(ey), std::move(ez), t, u, tMin, tMax, uMin, uMax);
            graph->setColor(getColor(color));
            emit graphChanged(this, graph);
        } else {
            qDebug() << "Unknown tab?";
        }
    } catch (const std::exception& e) {
        setErrorMsg(QLatin1String(e.what()));
    }
}

void Graph3DProperties::deletePressed() {
    deleteLater();
}

void Graph3DProperties::setErrorMsg(QString str) {
    errorMsg->setText(str);
}
