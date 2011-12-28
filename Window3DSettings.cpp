#include "Window3DSettings.h"

#include <QDebug>

Window3DSettings::Window3DSettings(QWidget* parent): QGroupBox(parent) {
    setupUi(this);
    xMin->setValidator(new QDoubleValidator());
    xMax->setValidator(new QDoubleValidator());
    yMin->setValidator(new QDoubleValidator());
    yMax->setValidator(new QDoubleValidator());
    zMin->setValidator(new QDoubleValidator());
    zMax->setValidator(new QDoubleValidator());
    lightX->setValidator(new QDoubleValidator());
    lightY->setValidator(new QDoubleValidator());
    lightZ->setValidator(new QDoubleValidator());
    connect(showAxesBox, SIGNAL(toggled(bool)), SIGNAL(showAxes(bool)));
}

void Window3DSettings::onChanged() {
    Vector3D boxa(xMin->text().toFloat(), yMin->text().toFloat(), zMin->text().toFloat());
    Vector3D boxb(xMax->text().toFloat(), yMax->text().toFloat(), zMax->text().toFloat());
    Vector3D lightPos(lightX->text().toFloat(), lightY->text().toFloat(), lightZ->text().toFloat());
    if (boxb.x() > boxa.x() && boxb.y() > boxa.y() && boxb.z() > boxa.z()) {
        emit boxChanged(boxa, boxb);
    }
    emit lightChanged(lightPos);
}
