#include "ParametricGraph3D.h"

#include <QtConcurrentRun>

#include <boost/config/suffix.hpp>

BOOST_CONSTEXPR_OR_CONST std::size_t ParametricGraph3D::numPts;

ParametricGraph3D::ParametricGraph3D(QObject* parent): Graph3D(parent),
        tPts(VECTOR_ALLOC(numPts)), uPts(VECTOR_ALLOC(numPts)) {
}

void ParametricGraph3D::reset(
        EPtr x, EPtr y, EPtr z,
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

void ParametricGraph3D::compute() {
    // Setup
    _mm_empty();
    
    Buffer3D buffer(width(), height(), transform());
    buffer.setColor(color().rgba());
    buffer.setLight(light());
    tDist = std::uniform_real_distribution<Number>(tMin, tMax);
    uDist = std::uniform_real_distribution<Number>(uMin, uMax);

    forever {
        // generate t,u
        if (cancelled()) return;
        for (int i = 0; i < numPts; ++i) {
            tPts[i] = tDist(engine);
            uPts[i] = uDist(engine);
        }
        
        // evaluate everything
        if (cancelled()) return;
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

        // calculate cross products to find normal vectors
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
            buffer.drawTransformLitPoint(Vector3D<float>(px[i], py[i], pz[i]), Vector3D<float>(ox[i], oy[i], oz[i]));
        }

        VECTOR_FREE(px);
        VECTOR_FREE(py);
        VECTOR_FREE(pz);
        VECTOR_FREE(ox);
        VECTOR_FREE(oy);
        VECTOR_FREE(oz);

        emit updated(new Buffer3D(buffer.copy()));
    }
}
