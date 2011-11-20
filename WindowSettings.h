#ifndef _WINDOWSETTINGS_H_
#define _WINDOWSETTINGS_H_

#include "ui_WindowSettings.h"
#include <QGroupBox>

class WindowSettings : public QGroupBox, private Ui_WindowSettings {
    Q_OBJECT
public:
    WindowSettings(QWidget* parent = 0);
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
