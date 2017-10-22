#ifndef VMID_SRC_EXHAUSTIVE_DECOMPOSER_HPP
#define VMID_SRC_EXHAUSTIVE_DECOMPOSER_HPP

#include <vector>
#include "decomposer.hpp"

class ExhaustiveDecomposer : public Decomposer
{
public:
    /**
     * @brief decompose()で使用するbaseとbasisを設定する。
     * @param base  2であれば二値分解、3であれば三値分解。
     * @param basis 基底行列の基底数。大きいほど近似精度が向上する。
     */
    explicit ExhaustiveDecomposer(const int base=2, const int basis=1);

    /**
     * @brief 重みベクトル/行列を近似するために整数基底分解を実行する。
     * @param[in]  weight    重みベクトル/行列
     * @param[out] basis_mat 二値/三値の基底行列
     * @param[out] coeff     スケール係数ベクトル/行列
     * @param[in]  n_inits   初期値変更回数
     */
    template <typename Type>
    void decompose(const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& weight,
                   Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic>      & basis_mat,
                   Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>      & coeff,
                   const int                                                  n_inits=20) const;

    template <typename Type>
    void decompose(const Eigen::Matrix<Type, Eigen::Dynamic, 1>             & weight,
                   Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic>      & basis_mat,
                   Eigen::Matrix<Type, Eigen::Dynamic, 1>                   & coeff,
                   const int                                                  n_inits=20) const;

    /**
     * @brief 重みベクトル/行列を近似するために整数基底分解と分解した基底行列を64bitの配列に圧縮する。
     * @param[in]  weight    重みベクトル/行列
     * @param[out] basis_mat 二値/三値の基底行列
     * @param[out] coeff     スケール係数ベクトル/行列
     * @param[in]  n_inits   初期値変更回数
     */
    template <typename Type>
    void decompose(const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& weight,
                   Eigen::Matrix<uint64_t, Eigen::Dynamic, Eigen::Dynamic>  & basis_mat,
                   Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>      & coeff,
                   const int                                                  n_inits=20) const;

    template <typename Type>
    void decompose(const Eigen::Matrix<Type, Eigen::Dynamic, 1>             & weight,
                   Eigen::Matrix<uint64_t, Eigen::Dynamic, Eigen::Dynamic>  & basis_mat,
                   Eigen::Matrix<Type, Eigen::Dynamic, 1>                   & coeff,
                   const int                                                  n_inits=20) const;

private:
    /**
     * @brief 基底行列とスケール係数ベクトル/行列を最適化する
     * @param[in]  weight    重み行列
     * @param[out] basis_mat 二値/三値の基底行列
     * @param[out] coeff     スケール係数行列x
     * @param[in]  candidates basis vector.
     */
    template <typename Type>
    void update(const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> & weight,
                Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic>       & basis_mat_opt,
                Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>       & coeff_opt,
                const std::vector<Eigen::Matrix<Type, 1, Eigen::Dynamic> >& candidates) const;
};

#include "exhaustive_decomposer_impl.hpp"

#endif
