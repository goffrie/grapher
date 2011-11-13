#include "GraphProperties.h"
#include "Parser.h"

GraphProperties::GraphProperties(Variable _x, Variable _y, QWidget* parent): QGroupBox(parent), x(_x), y(_y) {
    setupUi(this);
}

void GraphProperties::textChanged() {
    try {
        boost::unordered_map<std::string, Expression*> vars;
        vars.insert(std::make_pair<std::string, Expression*>("x", &x));
        vars.insert(std::make_pair<std::string, Expression*>("y", &y));
        Thing* t = Parser::parse(equation->text().toStdString(), vars);
        Equation* eqn = dynamic_cast<Equation*>(t);
        if (eqn == NULL) {
            delete t;
            setErrorMsg(QLatin1String("Didn't get an equation!"));
        } else {
            setErrorMsg(QString());
            emit equationChanged(this, eqn);
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
