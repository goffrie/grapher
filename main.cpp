#include <QApplication>

#include <iostream>

#include <gsl/gsl_errno.h>

#include "Expression.h"
#include "Grapher3D.h"
#include "Graph3D.h"
#include "MainWindow.h"
#include "Render3D.h"

#include "Parser.h"
#include "dynamic_unique_cast.h"

int main(int argc, char *argv[]) {
    gsl_set_error_handler_off();

    QApplication app(argc, argv);
    app.setApplicationName("Grapher");
    qRegisterMetaType<Variable>();
    qRegisterMetaType<Vector3D<float>>();
    qRegisterMetaType<Vector3D<Vc::float_v>>();
    qRegisterMetaType<Transform3D>();

    std::unordered_map<std::string, Expression*> vars;
    Variable x("x"), y("y"), z("z"), t("t");
    vars.insert(std::make_pair(std::string("x"), &x));
    vars.insert(std::make_pair(std::string("y"), &y));
    vars.insert(std::make_pair(std::string("z"), &z));
    vars.insert(std::make_pair(std::string("t"), &t));
    if (0) {
        std::string s = "x + y * z^4";
        x.id->type = y.id->type = z.id->type = Variable::Id::Constant;
        float* _x = new float(2), *_y = new float(3), *_z = new float(5);
        x.id->p = _x;
        y.id->p = _y;
        z.id->p = _z;
        std::unique_ptr<Thing> p = Parser::parse(s, vars);
        std::cout << p->toString() << std::endl;
        EPtr q = dynamic_unique_cast<Expression>(std::move(p));
        MathContext ctx = MathContext::defaultContext();
        WEvalFunc e = q->evaluator(ctx);
        double r = e();
        std::cout << r << std::endl;
        return 0;
    }
    while (false) {
        try {
            std::string s = "t * t ^ 3";
            //std::getline(std::cin, s);
            std::unique_ptr<Thing> p = Parser::parse(s, vars);
            std::cout << p->toString() << std::endl;
            EPtr q = dynamic_unique_cast<Expression>(std::move(p));
            if (q.get()) {
                q = q->expand();
                std::cout << q->toString() << std::endl;
                q = q->facsum(t);
                std::cout << q->toString() << std::endl;
            }
        } catch (std::exception& e) {
            std::cerr << e.what() << std::endl;
        }
        return 0;
    }

    MainWindow* window = new MainWindow();
    window->show();
#if 0
    Grapher3D* g = new Grapher3D();
//    g->setBox(Vector3D(-5,-2,-5),Vector3D(5,2,5));
//    g->setBox(Vector3D(-1,-1,-1),Vector3D(1,1,1));
    g->setBox(Vector3D(-2,-1,-2),Vector3D(2,1,2));
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
#endif
    return app.exec();
}
