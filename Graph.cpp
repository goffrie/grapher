#include "Graph.h"

#include <QDebug>

Graph::Graph(QObject* parent): QObject(parent), m_cancelled(0), m_thread(new GraphRunner(this)) {
}

Graph::~Graph() {
    Q_ASSERT(!m_thread->isRunning());
}

bool Graph::cancelled() const {
    return m_cancelled;
}

void Graph::restart() {
    Q_ASSERT(!m_thread->isRunning());
    m_cancelled = false;
    emit started();
    m_thread->start();
}

void Graph::stop() {
    m_cancelled = true;
    m_thread->wait();
    emit stopped();
}

void Graph::dispose() {
    m_cancelled = true;
    m_thread->wait();
    deleteLater();
}
