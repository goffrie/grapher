#ifndef _WINDOW3DSETTINGS_H_
#define _WINDOW3DSETTINGS_H_

#include "ui_Window3DSettings.h"

#include "Render3D.h"

class Window3DSettings : public QGroupBox, private Ui_Window3DSettings {
    Q_OBJECT
public:
    Window3DSettings(QWidget* parent = 0);
signals:
    void boxChanged(Vector3D<float> boxa, Vector3D<float> boxb);
    void lightChanged(Vector3D<float> lightPos);
    void showAxes(bool showAxes);
public slots:
    void onChanged();
    void setLight(Vector3D<float> lightPos);
};

#endif
