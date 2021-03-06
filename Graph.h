#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <QObject>
#include <QColor>
#include <QImage>
#include <QThread>

/// Abstract base class for graphs.
class Graph: public QObject {
    Q_OBJECT
    QColor m_color;
    bool m_cancelled;
    QThread* m_thread;
    friend class GraphRunner;
public:
    Graph(QObject* parent = 0);
    /// Destructor. Currently only performs a safety check.
    /// Do not destroy a Graph if its thread is running. If in doubt, call dispose() instead.
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
    /// Safely deletes the graph.
    void dispose();
signals:
    /// Emitted when the graph starts computing.
    void started();
    /// Emitted when the graph stops computing.
    void stopped();
};

/// Thread for running graphs.
class GraphRunner: public QThread {
    Q_OBJECT
    friend class Graph;
    Graph* graph;
    GraphRunner(Graph* g): QThread(g), graph(g) { }
    void run() override {
        graph->compute();
    }
};

#endif
