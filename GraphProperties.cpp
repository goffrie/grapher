#include "GraphProperties.h"
#include "Parser.h"
#include "dynamic_unique_cast.h"

GraphProperties::GraphProperties(Variable _x, Variable _y, QWidget* parent): QGroupBox(parent), x(_x), y(_y) {
    setupUi(this);
}

void GraphProperties::textChanged() {
    try {
        boost::unordered_map<std::string, Expression*> vars;
        vars.insert(std::make_pair<std::string, Expression*>("x", &x));
        vars.insert(std::make_pair<std::string, Expression*>("y", &y));
        auto eqn = dynamic_unique_cast<Equation>(Parser::parse(equation->text().toStdString(), vars));
        if (!eqn) {
            setErrorMsg(QLatin1String("Didn't get an equation!"));
        } else {
            setErrorMsg(QString());
            emit equationChanged(this, eqn.release());
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
