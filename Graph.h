#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <QObject>
#include <QColor>
#include <QImage>
#include <QFuture>

/// Abstract base class for graphs.
class Graph: public QObject {
    Q_OBJECT
    QColor m_color;
    bool m_cancelled;
    QFuture<void> m_future;
public:
    Graph(QObject* parent = 0);
    virtual ~Graph();

    /// Gets the graph's color.
    QColor color() const { return m_color; }
protected:
    /// Check if it has been requested that the computation be cancelled.
    bool cancelled() const;
    /// Calculate the graph.
    virtual void compute() = 0;
public slots:
    /// Starts the graph's computation from the beginning.
    void restart();
    /// Stop this graph's computation as soon as possible.
    /// Called when the graph is no longer needed.
    /// Waits for the thread to finish.
    void stop();
    /// Sets the graph's color.
    void setColor(QColor c) { m_color = c; }
signals:
    /// Emitted when the graph starts computing.
    void started();
    /// Emitted when the graph stops computing.
    void stopped();
};

#endif
