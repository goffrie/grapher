#include "Graph.h"

#include <QtConcurrentRun>

Graph::Graph(QObject* parent): QObject(parent), m_cancelled(0) {
}

Graph::~Graph() {
    m_cancelled = true;
    m_future.waitForFinished();
}

bool Graph::cancelled() const {
    return m_cancelled;
}

void Graph::restart() {
    m_cancelled = false;
    emit started();
    m_future = QtConcurrent::run(this, &Graph::compute);
}

void Graph::stop() {
    m_cancelled = true;
    m_future.waitForFinished();
    emit stopped();
}
