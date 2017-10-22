#ifndef VMID_SRC_DECOMPOSER_HPP
#define VMID_SRC_DECOMPOSER_HPP

#include <stdexcept>
#include "type.hpp"

/**
 * 整数基底分解を行うベースクラス。
 */
class Decomposer
{
public:
    Decomposer(const int base, const int basis)
    {
        set_base(base);
        set_basis(basis);
    }
    void set_base(const int base)
    {
        if (base != 2 && base != 3) {
            throw std::invalid_argument("``base`` must be 2 or 3.");
        }
        base_ = base;
    }
    void set_basis(const int basis)
    {
        if (basis <= 0) {
            throw std::invalid_argument("``basis`` must be 1 <= basis.");
        }
        basis_ = basis;
    }

protected:
    int base_;
    int basis_;
};

#endif
