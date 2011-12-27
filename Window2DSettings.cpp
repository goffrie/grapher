#include "Window2DSettings.h"

#include <QDebug>

Window2DSettings::Window2DSettings(QWidget* parent): QGroupBox(parent) {
    setupUi(this);
    xMin->setValidator(new QDoubleValidator());
    xMax->setValidator(new QDoubleValidator());
    yMin->setValidator(new QDoubleValidator());
    yMax->setValidator(new QDoubleValidator());
}

void Window2DSettings::onChanged() {
    qreal xmin = xMin->text().toDouble();
    qreal xmax = xMax->text().toDouble();
    qreal ymin = yMin->text().toDouble();
    qreal ymax = yMax->text().toDouble();
    bool good = true;
    if (xmax <= xmin) {
        
        good = false;
    }
    if (ymax <= ymin) {
        good = false;
    }
    if (good) {
        emit rectChanged(QRectF(xmin, ymin, xmax - xmin, ymax - ymin));
    }
}

void Window2DSettings::setShowAxes(bool _showAxes) {
    emit showAxes(_showAxes);
}

void Window2DSettings::setShowGrid(bool _showGrid) {
    emit showGrid(_showGrid);
}
