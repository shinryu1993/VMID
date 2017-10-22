#ifndef VMID_SRC_GREEDY_DECOMPOSER_IMPL_HPP
#define VMID_SRC_GREEDY_DECOMPOSER_IMPL_HPP

#include <cassert>
#include "algorithm.hpp"
#include "functor.hpp"
#include "plog/Log.h"

inline GreedyDecomposer::GreedyDecomposer(const int base, const int basis)
    : Decomposer(base, basis)
{
    if (base != 2) {
        throw std::invalid_argument("``base`` must be 2.");
    }
}

template <typename Type>
void GreedyDecomposer::decompose(
    const Eigen::Matrix<Type, Eigen::Dynamic, 1>       & weight,
    Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic>& basis_mat,
    Eigen::Matrix<Type, Eigen::Dynamic, 1>             & coeff)
{
    typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> VectorXx;
    assert(1 <= weight.rows());
    assert(1 == weight.cols());

    const int n_dims = weight.rows();
    basis_mat.resize(n_dims, basis_);
    coeff.resize(basis_);

    LOG_INFO << "Decomposing by greedy algorithm.";
    VectorXx residual = weight;
    for (int i = 0; i < basis_; ++i) {
        LOG_DEBUG << "basis number: " << i + 1;
        const VectorXx basis_vector = residual.unaryExpr(Sign<Type>());
        basis_mat.col(i) = basis_vector.template cast<int>();
        coeff(i) = (residual.transpose() * basis_vector)(0) / n_dims;
        residual -= basis_vector * coeff(i);
    }
    LOG_INFO << "Decomposed.";
}

template <typename Type>
void GreedyDecomposer::decompose(
    const Eigen::Matrix<Type, Eigen::Dynamic, 1>            & weight,
    Eigen::Matrix<uint64_t,  Eigen::Dynamic, Eigen::Dynamic>& basis_mat,
    Eigen::Matrix<Type, Eigen::Dynamic, 1>                  & coeff)
{
    Eigen::MatrixXi tmp_basis_mat;
    decompose<Type>(weight, tmp_basis_mat, coeff);
    compress(tmp_basis_mat, basis_mat);
}

#endif
