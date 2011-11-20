#ifndef _MAINWINDOW_H_
#define _MAINWINDOW_H_

#include "ui_MainWindow.h"

class MainWindow : public QMainWindow, private Ui_MainWindow {
    Q_OBJECT
public:
    MainWindow();
public slots:
    void newGraph();
    void scrollDown(QWidget* w);
};

#endif