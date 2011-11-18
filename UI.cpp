#include "UI.h"

#include <iostream>
#include "Parser.h"
#include "UIFactory.h"

UI::UI(UIFactory* factory, QWidget* parent): QSplitter(parent), m_factory(factory) {
    setupUi(this);

    {
        QSplitterHandle *handle = this->handle(1);
        QVBoxLayout *layout = new QVBoxLayout(handle);
        layout->setSpacing(0);
        layout->setMargin(0);

        QFrame *line = new QFrame(handle);
        line->setFrameShape(QFrame::VLine);
        line->setFrameShadow(QFrame::Sunken);
        layout->addWidget(line);
        layout->addWidget(line);
    }

    {
        verticalLayout->addWidget(m_factory->grapher());
        m_factory->grapher()->setMinimumWidth(150);
    }

    scrollAreaLayout->insertWidget(0, m_factory->windowSettings());

    fixPane();
}

void UI::newGraph() {
    QWidget* newGraph = m_factory->graphProperties();
    scrollAreaLayout->insertWidget(scrollAreaLayout->count()-2, newGraph);
    fixPane();
    QMetaObject::invokeMethod(this, "scrollDown", Qt::QueuedConnection, Q_ARG(QWidget*, newGraph)); // super hackiness
}

void UI::fixPane() {
    // <hackiness>
    scrollArea->takeWidget();
    scrollArea->setWidget(scrollAreaWidgetContents); // force size recalculation
    scrollAreaWidgetContents->setAutoFillBackground(false); // work around strangeness
    scrollAreaWidgetContents->updateGeometry();
    scrollArea->updateGeometry();
    updateGeometry();
    // </hackiness>
}


void UI::scrollDown(QWidget* w) {
    scrollArea->ensureWidgetVisible(w);
}
