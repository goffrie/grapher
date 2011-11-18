#ifndef _DYNAMIC_UNIQUE_CAST_H_
#define _DYNAMIC_UNIQUE_CAST_H_

template<typename B, typename A>
std::unique_ptr<B> dynamic_unique_cast(std::unique_ptr<A> a) {
    B* b = dynamic_cast<B*>(a.get());
    if (b != NULL) {
        a.release(); // let go of a's pointer; b holds it now
    } // otherwise, let a's destructor take care of destroying it
    return std::unique_ptr<B>(b);
}

template<typename B, typename A>
std::unique_ptr<B> dynamic_maybe_unique_cast(std::unique_ptr<A>& a) {
    B* b = dynamic_cast<B*>(a.get());
    if (b != NULL) {
        a.release(); // let go of a's pointer; b holds it now
    } // otherwise, let a keep holding on to it
    return std::unique_ptr<B>(b);
}

#endif