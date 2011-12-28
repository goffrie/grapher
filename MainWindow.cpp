#include "MainWindow.h"

#include "UI.h"
#include "UI2DFactory.h"
#include "UI3DFactory.h"

MainWindow::MainWindow() {
    setupUi(this);
    tabWidget->addTab(new UI(new UI2DFactory), "2D");
    tabWidget->addTab(new UI(new UI3DFactory), "3D");
    tabWidget->setCurrentIndex(0);
}
