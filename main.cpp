#include <QApplication>

#include <iostream>
#include "MainWindow.h"
#include "Grapher3D.h"
#include "Graph3D.h"
#include <gsl/gsl_errno.h>
#include "Expression.h"

int main(int argc, char *argv[]) {
    gsl_set_error_handler_off();
    
    QApplication app(argc, argv);
    app.setApplicationName("Grapher");
    qRegisterMetaType<Variable>();
    
    MainWindow* window = new MainWindow();
    window->show();
    Grapher3D* g = new Grapher3D();
//    g->setWindow(Vector3D(-5,-2,-5),Vector3D(5,2,5));
//    g->setWindow(Vector3D(-1,-1,-1),Vector3D(1,1,1));
    g->setWindow(Vector3D(-2,-1,-2),Vector3D(2,1,2));
    g->setLightSource(Vector3D(30,-30,30));
    QObject* o = new QObject;
    g->addGraph(o);
    ImplicitGraph3D* gr = new ImplicitGraph3D();
    Variable x(std::string("x")), y(std::string("y")), z(std::string("z"));
//    gr->reset(Equation::create(PowInt::create(x.ecopy(),2)+PowInt::create(y.ecopy(),2)+PowInt::create(z.ecopy(),2),Constant::create(1)), x, y, z);
//    gr->reset(Equation::create(y.ecopy(), Cos::create(x.ecopy()) + Sin::create(z.ecopy())), x, y, z);
//    gr->reset(Equation::create(y.ecopy(), Mul::create(Cos::create(x.ecopy()), Sin::create(z.ecopy()))), x, y, z);
    gr->reset(Equation::create(PowInt::create(PowInt::create(x.ecopy(), 2) + Constant::create(9/4.f) * PowInt::create(y.ecopy(), 2) + PowInt::create(z.ecopy(), 2) - Constant::create(1), 3) - PowInt::create(x.ecopy(), 2) * PowInt::create(-z.ecopy(), 3) - Constant::create(9/80.f) * PowInt::create(y.ecopy(), 2) * PowInt::create(-z.ecopy(), 3), Constant::create(0)), x, y, z);
    gr->setColor(QColor(230, 10, 10, 255));
    g->changeGraph(o, gr);
    g->show();
    return app.exec();
}
