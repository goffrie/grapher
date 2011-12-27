#include "MainWindow.h"

#include "UI.h"
#include "UI2DFactory.h"

MainWindow::MainWindow() {
    setupUi(this);
    tabWidget->addTab(new UI(new UI2DFactory), "2D");
}
