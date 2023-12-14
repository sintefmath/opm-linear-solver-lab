#define BOOST_TEST_MODULE TestReadBinary
#include <boost/test/unit_test.hpp>

#include <read_binary.hpp>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>


BOOST_AUTO_TEST_CASE(TestSPE1ReadMatrix)
{
    constexpr int blockSize = 3;
    using MatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double, blockSize, blockSize>>;

    auto matrixFromBinary = readBinaryAMGCLMatrix<MatrixType>("../examples/matrices/spe1/matrix.bin");

    MatrixType matrixFromMM;

    Dune::loadMatrixMarket(matrixFromMM, "../examples/matrices/spe1/matrix.mm");

    BOOST_REQUIRE_EQUAL(matrixFromBinary.N(), matrixFromMM.N());
    BOOST_REQUIRE_EQUAL(matrixFromBinary.nonzeroes(), matrixFromMM.nonzeroes());

    for (auto rowIterator = matrixFromMM.begin(); rowIterator != matrixFromMM.end(); ++rowIterator) {
        for (auto columnIterator = rowIterator->begin(); columnIterator != rowIterator->end(); ++columnIterator) {

            for (int i = 0; i < blockSize; ++i) {
                for (int j = 0; j < blockSize; ++j) {
                    BOOST_CHECK_EQUAL((*columnIterator)[i][j],
                                      matrixFromBinary[rowIterator.index()][columnIterator.index()][i][j]);
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(TestSPE1ReadVector)
{
    constexpr int blockSize = 3;
    using VectorType = Dune::BlockVector<Dune::FieldVector<double, blockSize>>;

    auto vectorFromBinary = readBinaryAMGCLVector<VectorType>("../examples/matrices/spe1/rhs.bin");

    VectorType vectorFromMM;

    Dune::loadMatrixMarket(vectorFromMM, "../examples/matrices/spe1/rhs.mm");

    BOOST_REQUIRE_EQUAL(vectorFromMM.N(), vectorFromBinary.N());

    for (size_t i = 0; i < vectorFromMM.N(); ++i) {
        for (size_t j = 0; j < blockSize; ++j) {
            BOOST_CHECK_EQUAL(vectorFromMM[i][j], vectorFromBinary[i][j]);
        }
    }
}
