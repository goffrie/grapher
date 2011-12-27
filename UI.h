#ifndef _UI_H_
#define _UI_H_

#include "ui_UI.h"

class UIFactory;

class UI : public QSplitter, private Ui_UI {
    Q_OBJECT
    UIFactory* m_factory;
public:
    explicit UI(UIFactory* factory, QWidget* parent = 0);
public slots:
    void newGraph();
    void fixPane();
    void scrollDown(QWidget* w);
};

#endif
