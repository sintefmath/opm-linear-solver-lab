#include <petscksp.h>
#include <petscconf.h>

#ifndef PETSC_HAVE_MPIUNI
#define HAVE_MPI 1
#endif

#include <boost/program_options.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>

Mat readMatrix(const std::string& path)
{
    Mat A;
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetFromOptions(A);

    Opm::Parallel::Communication comm(PETSC_COMM_WORLD);

    std::array<int,3> info;
    std::ifstream in;
    std::string line;

    OPM_BEGIN_PARALLEL_TRY_CATCH()
    if (comm.rank() == 0) {
        in.open(path);

        if (!in.good()) {
            throw std::runtime_error("Error loading matrix " + path);
        }

        std::getline(in, line);
        if (line != "%%MatrixMarket matrix coordinate real general") {
            throw std::runtime_error("Unexpected header in matrix " + path + ": " + line);
        }

        std::stringstream str;

        std::getline(in, line);
        str.str(line);
        str >> line >> line >> line >> info[0];

        if (line != "blocked") {
            throw std::runtime_error("Expect a blocked matrix");
        }

        std::getline(in, line);
        str.str(line);
        str >> info[1] >> info[1] >> info[2];
    }
    OPM_END_PARALLEL_TRY_CATCH("readMatrix", comm);

    comm.broadcast(info.data(), 3, 0);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, info[1], info[1]);
    MatMPIAIJSetPreallocation(A, info[2] / info[1], PETSC_NULLPTR, info[2] / info[1], PETSC_NULLPTR);
    MatSetBlockSizes(A, info[0], info[0]);
    MatSetUp(A);
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    if (comm.rank() == 0) {
        std::stringstream str;
        while (std::getline(in, line)) {
            str.clear();
            str.str(line);
            int i, j;
            double val;
            str >> i >> j >> val;
             MatSetValue(A, i-1, j-1, val, INSERT_VALUES);
        }
    }
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    return A;
}

Vec readVector(const std::string& path)
{
    Vec b;
    VecCreate(PETSC_COMM_WORLD, &b);
    VecSetFromOptions(b);

    Opm::Parallel::Communication comm(PETSC_COMM_WORLD);

    std::ifstream in(path);
    std::string line;
    std::array<int,2> info;

    OPM_BEGIN_PARALLEL_TRY_CATCH()
    if (comm.rank() == 0) {
        if (!in.good()) {
            throw std::runtime_error("Error loading vector " + path);
        }

        std::getline(in, line);
        if (line != "%%MatrixMarket matrix array real general") {
            throw std::runtime_error("Unexpected header in vector " + path + ": " + line);
        }

        std::stringstream str;
        int sanity;

        std::getline(in, line);
        str.str(line);
        str >> line >> line >> line >> info[0] >> sanity;

        if (line != "blocked" || sanity != 1) {
            throw std::runtime_error("Expect a blocked vector");
        }

        std::getline(in, line);
        str.clear();
        str.str(line);
        str >> info[1] >> sanity;
        if (sanity != 1) {
            throw std::runtime_error("Expect a vector");
        }
    }
    OPM_END_PARALLEL_TRY_CATCH("readVector", comm);

    comm.broadcast(info.data(), 2, 0);
    VecSetSizes(b, PETSC_DECIDE, info[1]);
    VecSetBlockSize(b, info[0]);
    VecSetFromOptions(b);
    VecSetUp(b);

    if (comm.rank() == 0) {
        int pos = 0;
        std::stringstream str;
        while (std::getline(in, line)) {
            str.clear();
            str.str(line);
            double val;
            str >> val;
            VecSetValue(b, pos++, val, INSERT_VALUES);
        }
    }

    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    return b;
}

int main(int argc, char** argv)
{
    namespace po = boost::program_options;
    po::options_description desc("Run matrix benchmark.");
    desc.add_options()("help", "Produce this help message.")(
        "matrix-file,m", po::value<std::string>()->required(), "Matrix filename.")(
        "initial-guess-file,x", po::value<std::string>()->default_value(""), "x (initial guess) filename.")(
        "rhs-file,y", po::value<std::string>()->required(), "y (right hand side) filename.");

    po::variables_map vm;

    try {
        po::store(po::command_line_parser(argc, argv).options(desc).allow_unregistered().run(), vm);

        po::notify(vm);
    } catch (const po::required_option& error) {
        std::cout << "Usage:\n\t" << argv[0]
                  << " -m <path to matrix file> -x <path to initial guess file> -y "
                     "<path to rhs file>"
                  << std::endl
                  << std::endl;

        std::cout << desc << std::endl;

        std::exit(EXIT_FAILURE);
    } catch (std::runtime_error& error) {
        std::cout << error.what() << std::endl;
        std::cout << "Usage:\n\t" << argv[0]
                  << " -m <path to matrix file> -x <path to initial guess file> -y "
                     "<path to rhs file>"
                  << std::endl
                  << std::endl;

        std::cout << desc << std::endl;

        std::exit(EXIT_FAILURE);
    }

    if (vm.count("help")) {
        std::cout << desc << "\n";
        std::exit(EXIT_FAILURE);
    }

    const auto matrixFilename = vm["matrix-file"].as<std::string>();
    const auto xFilename = vm["initial-guess-file"].as<std::string>();
    const auto rhsFilename = vm["rhs-file"].as<std::string>();

    PetscInitialize(&argc, &argv, 0, PETSC_NULLPTR);

    Mat A = readMatrix(matrixFilename);
    Vec b = readVector(rhsFilename);
    Vec x;
    if (!xFilename.empty())
        x = readVector(xFilename);
    else {
        VecDuplicate(b, &x);
    }

    KSP ksp;
    KSPCreate(PETSC_COMM_WORLD, &ksp);

    PC pc;
    KSPGetPC(ksp, &pc);
    PCSetFromOptions(pc);

    KSPSetOperators(ksp, A, A);
    if (!xFilename.empty()) {
        KSPSetInitialGuessNonzero(ksp, xFilename.empty() ? PETSC_FALSE : PETSC_TRUE);
    }

    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);

    Opm::Parallel::Communication comm(PETSC_COMM_WORLD);

    KSPSolve(ksp,b,x);
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp,&reason);
    if (reason < 0) {
        if (comm.rank() == 0) {
            std::cout << "Linear solve failed with reason " << KSPConvergedReasons[reason] << std::endl;
        }
        return 1;
    }

    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    if (comm.rank() == 0) {
        std::cout << "Success! Converged in " << its << " iterations" << std::endl;
    }

    MatDestroy(&A);
    VecDestroy(&b);
    VecDestroy(&x);
    KSPDestroy(&ksp);

    PetscFinalize();
    return 0;
}
