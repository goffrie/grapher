#include <QApplication>

#include <iostream>
#include "MainWindow.h"
#include <gsl/gsl_errno.h>

int main(int argc, char *argv[]) {
    gsl_set_error_handler_off();
	QApplication app(argc, argv);
	app.setApplicationName("Grapher");
    MainWindow* window = new MainWindow();
	window->show();
	return app.exec();
}
