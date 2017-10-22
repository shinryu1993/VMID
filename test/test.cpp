#include <cstdlib>
#include <ctime>
#include "plog/Log.h"
#include "plog/Appenders/ColorConsoleAppender.h"
#include "test_greedy_approximation.hpp"
#include "test_exhaustive_approximation.hpp"

void test_greedy_approximation(const Eigen::VectorXf fvector,
                               const Eigen::VectorXd dvector)
{
    {
        const int base = 2;
        const int basis = 5;

        LOG_INFO << "========================================";
        LOG_INFO << "Test of float type vector decompose by greedy algorithm.";
        greedy::test(base, basis, fvector);
        LOG_INFO << "========================================";

        LOG_INFO << "========================================";
        LOG_INFO << "Test of double type vector decompose by greedy algorithm.";
        greedy::test(base, basis, dvector);
        LOG_INFO << "========================================";
    }
}

void test_exhaustive_approximation(const Eigen::VectorXf& fvector,
                                   const Eigen::VectorXd& dvector,
                                   const Eigen::MatrixXf& fmatrix,
                                   const Eigen::MatrixXd& dmatrix)
{
    {
        const int base = 2;
        const int basis = 4;
        const int n_inits = 5;

        LOG_INFO << "========================================";
        LOG_INFO << "Test of float type vector decompose by exhaustive algorithm.";
        exhaustive::test(base, basis, n_inits, fvector);
        LOG_INFO << "========================================";
        LOG_INFO << "Test of double type vector decompose by exhaustive algorithm.";
        exhaustive::test(base, basis, n_inits, dvector);
        LOG_INFO << "========================================";
        LOG_INFO << "Test of float type matrix decompose by exhaustive algorithm.";
        exhaustive::test(base, basis, n_inits, fmatrix);
        LOG_INFO << "========================================";
        LOG_INFO << "Test of double type matrix decompose by exhaustive algorithm.";
        exhaustive::test(base, basis, n_inits, dmatrix);
        LOG_INFO << "========================================";
    }
    {
        const int base = 3;
        const int basis = 4;
        const int n_inits = 5;

        LOG_INFO << "========================================";
        LOG_INFO << "Test of float type vector decompose by exhaustive algorithm.";
        exhaustive::test(base, basis, n_inits, fvector);
        LOG_INFO << "========================================";
        LOG_INFO << "Test of double type vector decompose by exhaustive algorithm.";
        exhaustive::test(base, basis, n_inits, dvector);
        LOG_INFO << "========================================";
        LOG_INFO << "Test of float type matrix decompose by exhaustive algorithm.";
        exhaustive::test(base, basis, n_inits, fmatrix);
        LOG_INFO << "========================================";
        LOG_INFO << "Test of double type matrix decompose by exhaustive algorithm.";
        exhaustive::test(base, basis, n_inits, dmatrix);
        LOG_INFO << "========================================";
    }
}

int main(void)
{
    // 必須
    static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
    plog::init(plog::info, &consoleAppender);
    srand(static_cast<unsigned>(time(0)));

    const int n_dims = 10;
    const int n_outs = 5;
    const Eigen::VectorXf fvector = Eigen::VectorXf::Random(n_dims);
    const Eigen::VectorXd dvector = fvector.cast<double>();
    const Eigen::MatrixXf fmatrix = Eigen::MatrixXf::Random(n_dims, n_outs);
    const Eigen::MatrixXd dmatrix = fmatrix.cast<double>();
    test_greedy_approximation(fvector, dvector);
    test_exhaustive_approximation(fvector, dvector, fmatrix, dmatrix);

    return 0;
}
