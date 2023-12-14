#ifndef READ_BINARY_HPP
#define READ_BINARY_HPP
#include "config.hpp"
#include <fmt/format.h>
#include <string>

#if HAVE_AMGCL
#include <amgcl/io/binary.hpp>

template <class MatrixType>
MatrixType
readBinaryAMGCLMatrix(const std::string& filepath)
{
    std::vector<std::ptrdiff_t> ptr;
    std::vector<std::ptrdiff_t> col;
    std::vector<double> val;
    std::size_t nScalar = 0;

    amgcl::io::read_crs(filepath, nScalar, ptr, col, val);

    constexpr int blockRows = MatrixType::block_type::rows;
    constexpr int blockCols = MatrixType::block_type::cols;

    if (nScalar % blockRows != 0 || blockCols % blockCols != 0) {
        throw std::runtime_error(fmt::format(
            "Unexpected matrix size. Got n={}, blockRows = {}, blockCols = {}.", nScalar, blockRows, blockCols));
    }

    if (ptr.back() % (blockRows * blockCols) != 0) {
        throw std::runtime_error(fmt::format(
            "Unexpected nnz. Got ptr.back() = {}, blockRows = {}, blockCols = {}.", ptr.back(), blockRows, blockCols));
    }

    const auto n = nScalar / blockRows;
    const auto nnz = ptr.back() / (blockRows * blockCols);

    MatrixType matrix(n, n, nnz, MatrixType::row_wise);

    // Fill in sparsity pattern
    for (auto iter = matrix.createbegin(); iter != matrix.createend(); ++iter) {
        const auto scalarRowEnd = iter.index() * blockRows + blockRows;

        for (auto scalarRow = iter.index() * blockRows; scalarRow < scalarRowEnd; ++scalarRow) {
            for (auto column = ptr[scalarRow]; column < ptr[scalarRow + 1]; ++column) {
                const auto columnIndex = col[column];
                iter.insert(columnIndex / blockCols);
            }
        }
    }



    // Fill in values
    for (size_t row = 0; row < n; ++row) {
        for (size_t blockRow = 0; blockRow < blockRows; ++blockRow) {
            const auto scalarRow = row * blockRows + blockRow;

            for (auto scalarColumnPtr = ptr[scalarRow]; scalarColumnPtr < ptr[scalarRow + 1]; ++scalarColumnPtr) {
                const auto scalarColumn = col[scalarColumnPtr];
                const auto column = scalarColumn / blockCols;
                const auto blockColumn = scalarColumn % blockCols;
                matrix[row][column][blockRow][blockColumn] = val[scalarColumnPtr];
            }
        }
    }

    return matrix;
}

template <class VectorType>
VectorType
readBinaryAMGCLVector(const std::string& filepath)
{
    std::vector<double> val;
    std::size_t nScalar = 0;
    std::size_t mScalar = 0;

    amgcl::io::read_dense(filepath, nScalar, mScalar, val);
    constexpr int blockRows = VectorType::block_type::dimension;

    if (nScalar % blockRows != 0) {
        throw std::runtime_error(fmt::format("Unexpected vector size. Got n={}, blockRows = {}.", nScalar, blockRows));
    }

    if (mScalar != 1) {
        throw std::runtime_error(fmt::format("Got mScalar = {} (!= 1)", mScalar));
    }

    if (val.size() != nScalar) {
        throw std::runtime_error(
            fmt::format("Weird output from amgcl. Got val.size() = {}, but nScalar = {}.", val.size(), nScalar));
    }

    const auto n = nScalar / blockRows;

    VectorType vector(n);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < blockRows; ++j) {
            vector[i][j] = val[i * blockRows + j];
        }
    }

    return vector;
}

#else
template <class MatrixType>
MatrixType
readBinaryAMGCLMatrix(const std::string& filepath)
{
    throw std::runtime_error("You need to compile with amgcl to get binary support");
}

template <class VectorType>
VectorType
readBinaryAMGCLVector(const std::string& filepath)
{
    throw std::runtime_error("You need to compile with amgcl to get binary support");
}
#endif

#endif // READ_BINARY_HPP
