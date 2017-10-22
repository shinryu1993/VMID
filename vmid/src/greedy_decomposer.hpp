#ifndef VMID_SRC_GREEDY_DECOMPOSER_HPP
#define VMID_SRC_GREEDY_DECOMPOSER_HPP

#include "decomposer.hpp"

/**
 * Hare's proposed decomposition algorithm.
 * http://www.samhare.net/research/keypoints
 */
class GreedyDecomposer : public Decomposer
{
public:
    explicit GreedyDecomposer(const int base, const int basis);

    /**
     * @brief Integer decomposition for weight approximation.
     * @param[in]  weight    重みベクトル
     * @param[out] basis_mat 二値/三値基底行列
     * @param[out] coeff     スケール係数行列
     */
    template <typename Type>
    void decompose(const Eigen::Matrix<Type, Eigen::Dynamic, 1>       & weight,
                   Eigen::Matrix<int,  Eigen::Dynamic, Eigen::Dynamic>& basis_mat,
                   Eigen::Matrix<Type, Eigen::Dynamic, 1>             & coeff);

    /**
     * @brief Integer decomposition and compression for weight approximation.
     * @param[in]  weight    重みベクトル
     * @param[out] basis_mat 二値/三値基底行列。std::uint64_t
     * @param[out] coeff     スケール係数行列
     */
    template <typename Type>
    void decompose(const Eigen::Matrix<Type, Eigen::Dynamic, 1>            & weight,
                   Eigen::Matrix<uint64_t,  Eigen::Dynamic, Eigen::Dynamic>& basis_mat,
                   Eigen::Matrix<Type, Eigen::Dynamic, 1>                  & coeff);
};

#include "greedy_decomposer_impl.hpp"

#endif
