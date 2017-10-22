#ifndef VMID_SRC_FUNCTOR_HPP
#define VMID_SRC_FUNCTOR_HPP

#include <cstdlib>
#include "type.hpp"

/**
 * @brief サイン関数。0以上のとき1を返し、0未満のとき−1を返す。
 */
template <typename Type>
struct Sign {
    Type operator () (const Type v) const
    {
        if (static_cast<Type>(0) <= v) {
            return static_cast<Type>(1);
        } else {
            return static_cast<Type>(-1);
        }
    }
};

/**
 * @brief 2/3値で初期化する。
 */
struct BaseInitializer {
    explicit BaseInitializer(const int base)
        : base_(base)
    {}
    int operator () (const int v) const
    {
        return (Eigen::Vector3i() << -1, 1, 0).finished()(rand() % base_);
    }
    int base_;
};

#endif
