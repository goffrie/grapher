#include "MainWindow.h"

#include <iostream>
#include "Parser.h"
#include "GraphProperties.h"

MainWindow::MainWindow() {
    setupUi(this);
#if 0
    graph->setGeometry(100, 100, 200, 200);
    Variable x("x"), y("y");
    std::cout << x.toString() << ',' << y.toString() << std::endl;
    boost::unordered_map<std::string, Expression*> vars;
    vars.insert(std::make_pair<std::string, Expression*>("x", &x));
    vars.insert(std::make_pair<std::string, Expression*>("y", &y));
    Equation* eq = NULL;
    try {
        std::string str("x^2 + y^2 * (1 + 0.5 x) = 1");
        eq = dynamic_cast<Equation*>(Parser::parse(str, vars));
    } catch (std::exception& ex) {
        std::cout << ex.what() << std::endl;
        throw;
    }
//  Equation eq(new Add(new PowInt(x.copy(), 2), new PowInt(y.copy(), 2)), new Constant(4));
//  Equation eq(new PowInt(y.copy(), 2), new Add(new Sub(new PowInt(x.copy(), 3), x.copy()), new Constant(1)));
//  Equation eq(new PowInt(y.copy(), 6), new Sub(new PowInt(x.copy(), 2), new PowInt(x.copy(), 6)));
//  Equation eq(new Sub(new PowInt(y.copy(), 4), new PowInt(x.copy(), 4)), new Mul(x.copy(), y.copy()));
//  Equation eq(new Add(new PowInt(y.copy(), 3), new PowInt(x.copy(), 3)), new Mul(x.copy(), y.copy()));
/*    Equation eq(new PowInt(new Add(new PowInt(x.copy(), 2),new Sub(new Mul(x.copy(), y.copy()), new Constant(1))), 2),
                new Mul(new Sub(new Constant(1), new PowInt(x.copy(), 2)), new PowInt(new Sub(x.copy(), y.copy()), 2)));*/
/*    Equation eq(new Mul(new Add(new PowInt(x.copy(), 2), new PowInt(y.copy(), 2)), new Add(new PowInt(y.copy(), 2), new Add(x.copy(), new Constant(5)))),
                new Mul(new Constant(4*5), new Mul(x.copy(), new PowInt(y.copy(), 2))));*/
    
//  eq.a = new PowInt(Expression::Ptr(new Add(Expression::Ptr(new PowInt(x, 2)), Expression::Ptr(new PowInt(y, 2)))), 2);
//  eq.b = new Add(Expression::Ptr(new Mul(Expression::Ptr(new Constant(4)), Expression::Ptr(new PowInt(x, 2)))), Expression::Ptr(new Add(Expression::Ptr(new Mul(Expression::Ptr(new Constant(6)), Expression::Ptr(new PowInt(y, 2)))), Expression::Ptr(new Constant(2)))));
//  eq.a = new Mul(x, Expression::Ptr(new Add(Expression::Ptr(new PowInt(x, 2)), Expression::Ptr(new PowInt(y, 2)))));
//  eq.b = new Mul(Expression::Ptr(new Constant(4)), Expression::Ptr(new PowInt(y, 2)));
    graph->appendGraph(*eq, x, y);
    graph->resized();
    graph->setWindow(QRectF(-15, -15, 30, 30));
#endif
    QSplitterHandle *handle = splitter->handle(1);
    QVBoxLayout *layout = new QVBoxLayout(handle);
    layout->setSpacing(0);
    layout->setMargin(0);

    QFrame *line = new QFrame(handle);
    line->setFrameShape(QFrame::VLine);
    line->setFrameShadow(QFrame::Sunken);
    layout->addWidget(line);

    windowSettings->onChanged();
}

void MainWindow::newGraph() {
    GraphProperties* newGraph = new GraphProperties();
    graphsLayout->addWidget(newGraph);
    graph->addGraph(newGraph);
    connect(newGraph, SIGNAL(graphChanged(QObject*, Graph*)), graph, SLOT(changeGraph(QObject*, Graph*)));
}

