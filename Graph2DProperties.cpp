#include "Graph2DProperties.h"

#include "Parser.h"
#include "Expression.h"
#include "Graph2D.h"

#include <QDebug>
#include <exception>
#include <set>
#include <algorithm>

#include "dynamic_unique_cast.h"

int nextColor = 0;
std::set<int> recoveredColors;

QColor getColor(int c) {
    const static qreal factor = -(3 - std::sqrt(5)) / 2;
    const qreal fc = factor * c;
    return QColor::fromHsvF(fc - std::floor(fc), 0.9, 0.8);
}

Graph2DProperties::Graph2DProperties(QWidget* parent): QGroupBox(parent) {
    setupUi(this);
    par_tMin->setValidator(new QDoubleValidator());
    par_tMax->setValidator(new QDoubleValidator());
    if (recoveredColors.size() > 0) {
        color = *recoveredColors.begin();
        recoveredColors.erase(recoveredColors.begin());
    } else {
        color = nextColor++;
    }
    QPalette p = colorButton->palette();
    p.setColor(QPalette::Button, getColor(color));
    colorButton->setPalette(p);
}

Graph2DProperties::~Graph2DProperties() {
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

ParametricGraph* parametrize(std::unique_ptr<Equation>& eqn, Variable x, Variable y) {
    std::set<Variable> left, right;
    eqn->a->variables(left);
    eqn->b->variables(right);
/*    std::vector<Variable> intersection;
    std::set_intersection(left.begin(), left.end(), right.begin(), right.end(), std::back_insert_iterator(intersection));
    if (intersection.size() != 0) {
        // don't know how to deal with that
        return NULL;
    }*/
    Variable* l = dynamic_cast<Variable*>(eqn->a.get());
    Variable* r = dynamic_cast<Variable*>(eqn->b.get());
    ParametricGraph* ret = NULL;
    if (l != NULL && right.find(*l) == right.end()) {
        if (*l == x) {
            ret = new ParametricGraph;
            ret->reset(std::move(eqn->b), Variable::create(y), y, 0, 0);
            eqn.reset();
        } else if (*l == y) {
            ret = new ParametricGraph;
            ret->reset(Variable::create(x), std::move(eqn->b), x, 0, 0);
            eqn.reset();
        }
    } else if (r != NULL && left.find(*r) == left.end()) {
        if (*r == x) {
            ret = new ParametricGraph;
            ret->reset(std::move(eqn->a), Variable::create(y), y, 0, 0);
            eqn.reset();
        } else if (*r == y) {
            ret = new ParametricGraph;
            ret->reset(Variable::create(x), std::move(eqn->a), x, 0, 0);
            eqn.reset();
        }
    }
    return ret;
}

void Graph2DProperties::textChanged() {
    try {
        std::unordered_map<std::string, Expression*> vars;
        Constant _pi(M_PI);
        Constant _e(M_E);
        vars.insert(std::make_pair("pi", &_pi));
        vars.insert(std::make_pair("e", &_e));
        QWidget* current = graphTypeSwitcher->currentWidget();
        if (current == relation_tab) {
            Variable x("x"), y("y");
            vars.insert(std::make_pair<std::string, Expression*>("x", &x));
            vars.insert(std::make_pair<std::string, Expression*>("y", &y));
            std::unique_ptr<Thing> thing(Parser::parse(rel_equation->text().toStdString(), vars));
            std::unique_ptr<Equation> eqn = dynamic_maybe_unique_cast<Equation>(thing);
            if (eqn) {
                setErrorMsg(QString());
                ParametricGraph* para = parametrize(eqn, x, y);
                if (para) {
                    para->setColor(getColor(color));
                    emit graphChanged(this, para);
                } else {
                    ImplicitGraph* graph = new ImplicitGraph;
                    graph->reset(std::move(eqn), x, y);
                    graph->setColor(getColor(color));
                    emit graphChanged(this, graph);
                }
            } else {
                auto ineq = dynamic_unique_cast<Inequality>(std::move(thing));
                if (!ineq) {
                    throw InvalidInputException("Didn't get an equation!");
                }
                setErrorMsg(QString());
                InequalityGraph* graph = new InequalityGraph;
                graph->reset(std::move(ineq), x, y);
                graph->setColor(getColor(color));
                emit graphChanged(this, graph);
            }
        } else if (current == parametric_tab) {
            double tMin = par_tMin->text().toDouble(), tMax = par_tMax->text().toDouble();
            if (tMax <= tMin || !par_tMin->text().size() || !par_tMax->text().size()) {
                throw InvalidInputException("Bad interval!");
            }
            Variable t("t");
            vars.insert(std::make_pair<std::string, Expression*>("t", &t));
            auto ex = dynamic_unique_cast<Expression>(Parser::parse(par_x->text().toStdString(), vars));
            if (!ex) {
                throw InvalidInputException("x(t) isn't a function!");
            }
            auto ey = dynamic_unique_cast<Expression>(Parser::parse(par_y->text().toStdString(), vars));
            if (!ey) {
                throw InvalidInputException("y(t) isn't a function!");
            }
            setErrorMsg(QString());
            ParametricGraph* graph = new ParametricGraph;
            graph->reset(std::move(ex), std::move(ey), t, tMin, tMax);
            graph->setColor(getColor(color));
            emit graphChanged(this, graph);
        } else {
            qDebug() << "Unknown tab?";
        }
    } catch (const std::exception& e) {
        setErrorMsg(QLatin1String(e.what()));
    }
}

void Graph2DProperties::deletePressed() {
    deleteLater();
}

void Graph2DProperties::setErrorMsg(QString str) {
    errorMsg->setText(str);
}
