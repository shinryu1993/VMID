#ifndef VMID_TEST_TEST_EXHAUSTIVE_APPROXIMATION_HPP
#define VMID_TEST_TEST_EXHAUSTIVE_APPROXIMATION_HPP

#include "src/exhaustive_decomposer.hpp"

namespace exhaustive {

template <typename Type>
void test(const int base, const int limit_basis, const int n_inits,
          const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& subject)
{
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXx;
    for (int basis = 1; basis <= limit_basis; ++basis) {
        Eigen::MatrixXi basis_mat;
        MatrixXx coeff;
        ExhaustiveDecomposer decomposer(base, basis);
        decomposer.decompose(subject, basis_mat, coeff, n_inits);
        const Type squared_error = (subject - basis_mat.cast<Type>() * coeff).squaredNorm();
        LOG_INFO << "base: " << base << ", basis: " << basis << ", squared_error: " << squared_error;
    }
}

template <typename Type>
void test(const int base, const int limit_basis, const int n_inits,
          const Eigen::Matrix<Type, Eigen::Dynamic, 1>& subject)
{
    typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXx;
    MatrixXx subject_(Eigen::Map<const MatrixXx>(subject.data(), subject.rows(), subject.cols()));
    test(base, limit_basis, n_inits, subject_);
}

}  // namespace exhaustive

#endif
