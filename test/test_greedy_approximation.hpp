#ifndef VMID_TEST_TEST_GREEDY_APPROXIMATION_HPP
#define VMID_TEST_TEST_GREEDY_APPROXIMATION_HPP

#include "src/greedy_decomposer.hpp"

namespace greedy {

template <typename Type>
void test(const int base, const int limit_basis,
          const Eigen::Matrix<Type, Eigen::Dynamic, 1>& subject)
{
    typedef Eigen::Matrix<Type, Eigen::Dynamic, 1> VectorXx;
    for (int basis = 1; basis <= limit_basis; ++basis) {
        Eigen::MatrixXi basis_mat;
        VectorXx coeff;
        GreedyDecomposer decomposer(base, basis);
        decomposer.decompose(subject, basis_mat, coeff);
        const Type squared_error = (subject - basis_mat.cast<Type>() * coeff).squaredNorm();
        LOG_INFO << "base: " << base << ", basis: " << basis << ", squared_error: " << squared_error;
    }
}

}  // namespace greedy

#endif
