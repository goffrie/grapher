#include "Graph3D.h"

#include <QtConcurrentRun>
#include <QtConcurrentMap>

#include <boost/config/suffix.hpp>

Graph3D::Graph3D(QObject* parent) : Graph(parent) {

}

void Graph3D::setupRestart(const Transform3D& t,
        int width, int height,
        Vector3D<float> boxa,
        Vector3D<float> boxb,
        Vector3D<float> light) {
    cancel();
    m_width = width;
    m_height = height;
    m_a->boxa = boxa;
    m_a->boxb = boxb;
    m_a->light = light;
    m_a->transform = t;
    m_buf = Buffer3D(m_width, m_height, m_a->transform);
    findEyeRay();
    if (m_width > 0 && m_height > 0) {
        startThread();
    }
}

void Graph3D::findEyeRay() {
    Transform3D inv = m_a->transform.inverted();
    Vector3D<float> pt1(0.f,0.f,1.f);
    Vector3D<float> pt2(0.f,0.f,0.f);
    pt1 = inv * pt1;
    pt2 = inv * pt2;
    m_a->eyeray = (pt2 - pt1).normalized<6>();
}


BOOST_CONSTEXPR_OR_CONST std::size_t ParametricGraph3D::numPts;

ParametricGraph3D::ParametricGraph3D(QObject* parent): Graph3D(parent),
		tPts(VECTOR_ALLOC(numPts)), uPts(VECTOR_ALLOC(numPts)) {
}

void ParametricGraph3D::reset(
		std::unique_ptr<Expression> x,
		std::unique_ptr<Expression> y,
		std::unique_ptr<Expression> z,
		const Variable& t, const Variable& u,
		Number tMin, Number tMax,
		Number uMin, Number uMax) {
	this->x = x->simplify();
	this->y = y->simplify();
	this->z = z->simplify();
	x_t = x->derivative(t)->simplify();
	y_t = y->derivative(t)->simplify();
	z_t = z->derivative(t)->simplify();
	x_u = x->derivative(u)->simplify();
	y_u = y->derivative(u)->simplify();
	z_u = z->derivative(u)->simplify();
	this->t = t;
	this->u = u;
	this->t.id->type = Variable::Id::Vector;
	this->t.id->p = tPts.get();
	this->u.id->type = Variable::Id::Vector;
	this->u.id->p = uPts.get();
    this->tMin = tMin;
	this->tMax = tMax;
	this->uMin = uMin;
	this->uMax = uMax;
}

void ParametricGraph3D::cancel() {
    cancelled = true;
    future.waitForFinished();
    cancelled = false;
}

void ParametricGraph3D::startThread() {
    future = QtConcurrent::run(this, &ParametricGraph3D::restart);
}

void ParametricGraph3D::restart() {
    // setup
    _mm_empty();
    m_buf.setColor(m_color.rgba());
    m_buf.setLight(m_a->light);
    tDist = std::uniform_real_distribution<Number>(tMin, tMax);
    uDist = std::uniform_real_distribution<Number>(uMin, uMax);

begin:

    // generate t,u
    if (cancelled) return;
    for (int i = 0; i < numPts; ++i) {
        tPts[i] = tDist(engine);
        uPts[i] = uDist(engine);
    }
    
    // evaluate everything
    if (cancelled) return;
    QFuture<Vector> fx = QtConcurrent::run(x.get(), &Expression::evaluateVector, numPts);
    QFuture<Vector> fy = QtConcurrent::run(y.get(), &Expression::evaluateVector, numPts);
    QFuture<Vector> fz = QtConcurrent::run(z.get(), &Expression::evaluateVector, numPts);
    QFuture<Vector> fxt = QtConcurrent::run(x_t.get(), &Expression::evaluateVector, numPts);
    QFuture<Vector> fyt = QtConcurrent::run(y_t.get(), &Expression::evaluateVector, numPts);
    QFuture<Vector> fzt = QtConcurrent::run(z_t.get(), &Expression::evaluateVector, numPts);
    QFuture<Vector> fxu = QtConcurrent::run(x_u.get(), &Expression::evaluateVector, numPts);
    QFuture<Vector> fyu = QtConcurrent::run(y_u.get(), &Expression::evaluateVector, numPts);
    QFuture<Vector> fzu = QtConcurrent::run(z_u.get(), &Expression::evaluateVector, numPts);
    fx.waitForFinished();
    fy.waitForFinished();
    fz.waitForFinished();
    fxt.waitForFinished();
    fyt.waitForFinished();
    fzt.waitForFinished();
    fxu.waitForFinished();
    fyu.waitForFinished();
    fzu.waitForFinished();

    // calculate cross products
    VectorR pxt = fxt.result(), pyt = fyt.result(), pzt = fzt.result(),
            pxu = fxu.result(), pyu = fyu.result(), pzu = fzu.result(),
            ox = VECTOR_ALLOC(numPts),
            oy = VECTOR_ALLOC(numPts),
            oz = VECTOR_ALLOC(numPts);
#define size numPts
    VECTOR_LOOP {
        V(ox) = V(pyt) * V(pzu) - V(pzt) * V(pyu);
        V(oy) = V(pzt) * V(pxu) - V(pxt) * V(pzu);
        V(oz) = V(pxt) * V(pyu) - V(pyt) * V(pxu);
    }
#undef size
    VECTOR_FREE(pxt);
    VECTOR_FREE(pyt);
    VECTOR_FREE(pzt);
    VECTOR_FREE(pxu);
    VECTOR_FREE(pyu);
    VECTOR_FREE(pzu);

    // draw points
    VectorR px = fx.result(), py = fy.result(), pz = fz.result();
    for (std::size_t i = 0; i < numPts; ++i) {
        m_buf.drawTransformLitPoint(Vector3D<float>(px[i], py[i], pz[i]), Vector3D<float>(ox[i], oy[i], oz[i]));
    }

    VECTOR_FREE(px);
    VECTOR_FREE(py);
    VECTOR_FREE(pz);
    VECTOR_FREE(ox);
    VECTOR_FREE(oy);
    VECTOR_FREE(oz);

    emit updated();
    goto begin;
}
