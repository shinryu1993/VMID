#ifndef VMID_SRC_ALGORITHM_HPP
#define VMID_SRC_ALGORITHM_HPP

#include <cmath>
#include <vector>
#include "type.hpp"
#include "plog/Log.h"

/**
 * @brief 基底行列を最適化する際の候補ベクトルの生成。
 * @par 基底行列を全パターン作成するとメモリが足りなくなる可能性が高い。
 *      各行毎に必要な候補は同じであるため、1行分の候補を作成し各要素の計算時に使い回す。
 * @param[in]  base       分解するときの基底値。
 * @param[in]  basis      基底行列の基底数。
 * @param[out] candidates 生成した候補ベクトルを格納する。ただし、列ベクトルではなく行ベクトル。
 */
template <typename Type>
void generate_candidates(const int base, const int basis,
                         std::vector<Eigen::Matrix<Type, 1, Eigen::Dynamic> >& candidates)
{
    LOG_INFO << "Generating candidates of basis matrix.";
    typedef Eigen::Matrix<Type, 1, Eigen::Dynamic> RowVectorXx;
    const int n_candidates = static_cast<int>(std::pow(base, basis));
    candidates.resize(n_candidates);
    for (int i = 0; i < n_candidates; ++i) {
        // Decimal to binary digits.
        int value = i;
        RowVectorXx candidate(basis);
        for (int k = 1; k <= basis; ++k) {
            candidate(basis - k) = static_cast<Type>(value % base);
            value /= base;
        }
        if (base == 3) {
            candidates[i] = candidate.array() * 1.0 - 1.0;  // Convert {0, 1, 2} to {-1, 0, 1}
        } else {
            candidates[i] = candidate.array() * 2.0 - 1.0;  // Convert {0, 1}    to {-1, 1}
        }
    }
    LOG_INFO << "Generated candidates.";
    LOG_INFO << "Number of generated candidates: " << candidates.size();
}

/**
 * @brief 分解した基底行列を基底毎に64ビット変数(unsigned long long)に圧縮する(詰め込む)。
 * @param[in]  src 分解された圧縮前の基底行列。
 * @param[out] dst 圧縮したbasis_matを格納する基底行列。
 */
inline void compress(const Eigen::MatrixXi& src, MatrixXl& dst)
{
    LOG_INFO << "Compressing basis matrix.";
    static const int VAR_BITS = sizeof(uint64_t) * 8;  // byte to bit
    const int n_dims = src.rows();
    const int basis = src.cols();
    const int n_compressed_dims = std::ceil(static_cast<float>(n_dims) / VAR_BITS);
    dst = MatrixXl::Zero(n_compressed_dims, basis);

    // compress
    for (int basis_no = 0; basis_no < basis; ++basis_no) {
        for (int dim_no = 0; dim_no < n_dims; ++dim_no) {
            const int index = dim_no / VAR_BITS;
            dst(index, basis_no) <<= 1;
            if (src(dim_no, basis_no) == 1) {
                dst(index, basis_no) += 1;
            }
        }
    }
    LOG_INFO << "Compressed basis matrix.";
}

#endif
