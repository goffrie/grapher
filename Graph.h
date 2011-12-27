#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <QObject>
#include <QTransform>
#include <QColor>

class QImage;

class Graph : public QObject {
    Q_OBJECT
public:
    Graph(QObject* parent = 0);
    virtual ~Graph() { }

    virtual void cancel() = 0;
    virtual void setColor(QColor c) { m_color = c; }
protected:
    virtual void startThread() = 0;
    int m_width, m_height;
    QColor m_color;
signals:
    void updated();
};

#endif
