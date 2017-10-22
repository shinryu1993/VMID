#ifndef VMID_SRC_EXHAUSTIVE_DECOMPOSER_IMPL_HPP
#define VMID_SRC_EXHAUSTIVE_DECOMPOSER_IMPL_HPP

#include <cassert>
#include <limits>
#include <vector>
#include "algorithm.hpp"
#include "functor.hpp"
#include "eigen3/Eigen/LU"
#include "plog/Log.h"

inline ExhaustiveDecomposer::ExhaustiveDecomposer(const int base, const int basis)
    : Decomposer(base, basis)
{}

/**
 * Weight approximation.
 */
template <typename Type>
void ExhaustiveDecomposer::decompose(
    const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& weight,
    Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic>      & basis_mat_opt,
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>      & coeff_opt,
    const int                                                  n_inits) const
{
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXx;
    typedef Eigen::Matrix<Type, 1, Eigen::Dynamic> RowVectorXx;
    using Eigen::MatrixXi;

    assert(0 < n_inits);
    assert(0 < weight.size());

    // 基底行列の候補を作成
    std::vector<RowVectorXx> candidates;
    generate_candidates(base_, basis_, candidates);

    // 分解開始
    LOG_INFO << "Decomposing by exhaustive algorithm.";
    const int n_dims = weight.rows();
    Type min_cost = std::numeric_limits<Type>::max();
    for (int i = 0; i < n_inits; ++i) {
        MatrixXi basis_mat = MatrixXi(n_dims, basis_).unaryExpr(BaseInitializer(base_));
        MatrixXx coeff;
        update(weight, basis_mat, coeff, candidates);

        const MatrixXx approximated = basis_mat.cast<Type>() * coeff;
        const Type cost = (weight - approximated).squaredNorm();  // Vector: L2 Norm, Matrix: Frobenius Norm
        if (cost < min_cost) {
            min_cost = cost;
            basis_mat_opt = basis_mat;
            coeff_opt = coeff;
        }
        LOG_DEBUG << "Init number: " << i + 1 << " / " << n_inits
                  << "MSE: " << cost << std::endl;
    }
    LOG_INFO << "Decomposed.";
}

template <typename Type>
void ExhaustiveDecomposer::decompose(
    const Eigen::Matrix<Type, Eigen::Dynamic, 1>       & weight,
    Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic>& basis_mat,
    Eigen::Matrix<Type, Eigen::Dynamic, 1>             & coeff,
    const int                                            n_inits) const
{
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXx;
    typedef Eigen::Matrix<Type, Eigen::Dynamic, 1>              VectorXx;
    MatrixXx weight_(Eigen::Map<const MatrixXx>(weight.data(), weight.rows(), weight.cols()));
    MatrixXx coeff_(Eigen::Map<MatrixXx>(coeff.data(), coeff.rows(), coeff.cols()));
    decompose(weight_, basis_mat, coeff_, n_inits);
    coeff = VectorXx(coeff_);
}

/**
 * Weight approximation and compression.
 */
template <typename Type>
void ExhaustiveDecomposer::decompose(
    const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& weight,
    Eigen::Matrix<uint64_t,  Eigen::Dynamic, Eigen::Dynamic> & basis_mat,
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>      & coeff,
    const int                                                  n_inits) const
{
    Eigen::MatrixXi tmp_basis_mat;
    decompose<Type>(weight, tmp_basis_mat, coeff, n_inits);
    compress(tmp_basis_mat, basis_mat);
}

template <typename Type>
void ExhaustiveDecomposer::decompose(
    const Eigen::Matrix<Type, Eigen::Dynamic, 1>            & weight,
    Eigen::Matrix<uint64_t,  Eigen::Dynamic, Eigen::Dynamic>& basis_mat,
    Eigen::Matrix<Type, Eigen::Dynamic, 1>                  & coeff,
    const int                                                 n_inits) const
{
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXx;
    typedef Eigen::Matrix<Type, Eigen::Dynamic, 1>              VectorXx;
    using Eigen::MatrixXi;
    MatrixXx weight_(Eigen::Map<const MatrixXx>(weight.data(), weight.rows(), weight.cols()));
    MatrixXx coeff_(Eigen::Map<MatrixXx>(coeff.data(), coeff.rows(), coeff.cols()));
    MatrixXi tmp_basis_mat;
    decompose<Type>(weight_, tmp_basis_mat, coeff_, n_inits);
    compress(tmp_basis_mat, basis_mat);
    coeff = VectorXx(coeff_);
}

/**
 * Updating basis matrix and vector/matrix coefficient.
 */
template <typename Type>
void ExhaustiveDecomposer::update(
    const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> & weight,
    Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic>       & basis_mat_opt,
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>       & coeff_opt,
    const std::vector<Eigen::Matrix<Type, 1, Eigen::Dynamic> >& candidates) const
{
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXx;
    typedef Eigen::Matrix<Type, 1, Eigen::Dynamic>              RowVectorXx;
    int converged_counter = 0;
    Type min_cost = std::numeric_limits<Type>::max();

    MatrixXx basis_mat = basis_mat_opt.cast<Type>();
    while (true) {
        // 最小二乗法によりスケール係数ベクトルを最適化する
        const MatrixXx basis_mat_t = basis_mat.transpose();
        const MatrixXx coeff = (basis_mat_t * basis_mat).fullPivLu().solve(basis_mat_t * weight);
        assert(coeff.rows() == basis_);
        assert(coeff.cols() == weight.cols());

        // 基底行列を網羅的に探索し最適化する
        const int n_candidates = candidates.size();
        std::vector<RowVectorXx> approximated(n_candidates);
        for (int i = 0; i < n_candidates; ++i) {
            approximated[i] = candidates[i] * coeff;
        }
        const int n_dims = weight.rows();
        for (int dim_no = 0; dim_no < n_dims; ++dim_no) {
            Type candidates_min_error = std::numeric_limits<Type>::max();
            for (int cand_no = 0; cand_no < n_candidates; ++cand_no) {
                const Type error = (weight.row(dim_no) - approximated[cand_no]).squaredNorm();
                if (error < candidates_min_error) {
                    candidates_min_error = error;
                    basis_mat.row(dim_no) = candidates[cand_no];
                }
            }
        }

        // 誤差が収束していればカウントする
        assert(weight.rows() == basis_mat.rows());
        assert(weight.cols() == coeff.cols());
        assert(basis_mat.cols() == coeff.rows());
        const Type cost = (weight - basis_mat * coeff).squaredNorm();
        if (std::fabs(cost - min_cost) <= 1e-10) {
            ++converged_counter;
        } else {
            converged_counter =  0;
        }

        if (cost < min_cost) {
            min_cost = cost;
            basis_mat_opt = basis_mat.template cast<int>();
            coeff_opt = coeff;
        }
        if (10 <= converged_counter || min_cost <= 1e-10) return;
        LOG_DEBUG << "cost: " << std::setprecision(std::numeric_limits<Type>::digits10 + 1) << cost;
    }
}

#endif
