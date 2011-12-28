#ifndef _WINDOW2DSETTINGS_H_
#define _WINDOW2DSETTINGS_H_

#include "ui_Window2DSettings.h"
#include <QGroupBox>

class Window2DSettings : public QGroupBox, private Ui_Window2DSettings {
    Q_OBJECT
public:
    Window2DSettings(QWidget* parent = 0);
signals:
    void rectChanged(QRectF rect);
    void showAxes(bool showAxes);
    void showGrid(bool showGrid);
public slots:
    void setShowAxes(bool showAxes);
    void setShowGrid(bool showGrid);
    void onChanged();
};

#endif
