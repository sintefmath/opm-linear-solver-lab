#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#pragma GCC push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <opm/simulators/linalg/cuistl/CuSeqILU0.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/PreconditionerAdapter.hpp>
#pragma GCC pop
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <opm/simulators/linalg/ilufirstelement.hh>
#include <opm/simulators/linalg/matrixblock.hh>

#include <fmt/format.h>

#include <limits>
#include <memory>
#include <random>

#include <boost/program_options.hpp>

template <class X>
class OnlyPreconditionSolver : public Dune::IterativeSolver<X, X>
{
public:
    using typename Dune::IterativeSolver<X, X>::domain_type;
    using typename Dune::IterativeSolver<X, X>::range_type;
    using typename Dune::IterativeSolver<X, X>::field_type;
    using typename Dune::IterativeSolver<X, X>::real_type;

    void apply(X& x, X& b, double reduction, Dune::InverseOperatorResult& res) override
    {
        return apply(x, b, res);
    }

    // copy base class constructors
    using Dune::IterativeSolver<X, X>::IterativeSolver;
    void apply(X& x, X& b, Dune::InverseOperatorResult& res) override
    {
        _prec->pre(x, b); // prepare preconditioner
        _prec->apply(x, b);
        _prec->post(x);

        res.converged = true;
        res.iterations = 1;
    }

protected:
    using Dune::IterativeSolver<X, X>::_prec;
};

template <int dim, template <class> class Solver = Dune::BiCGSTABSolver, class T = double>
std::tuple<unsigned long long, Dune::InverseOperatorResult, bool>
readAndSolveCPU(const auto jsonConfigCPUFilename,
                const auto xFilename,
                const auto matrixFilename,
                const auto rhsFilename)
{
    using M = Opm::MatrixBlock<T, dim, dim>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<T, dim>>;
    using CuILU0 = Opm::cuistl::CuSeqILU0<SpMatrix, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>;
    using Operator = Dune::MatrixAdapter<SpMatrix, Vector, Vector>;
    using PrecFactory = Opm::PreconditionerFactory<Operator, Dune::Amg::SequentialInformation>;

    Opm::PropertyTree configurationCPU(jsonConfigCPUFilename);
    bool transpose = false;
    if (configurationCPU.get<std::string>("preconditioner.type") == "cprt") {
        transpose = true;
    }

    SpMatrix B;
    Vector x, rhs;

    using MatrixDuneIsExpecting = Dune::BCRSMatrix<Dune::FieldMatrix<double, dim, dim>>;

    Dune::loadMatrixMarket(reinterpret_cast<MatrixDuneIsExpecting&>(B), matrixFilename);
    Dune::loadMatrixMarket(x, xFilename);
    Dune::loadMatrixMarket(rhs, rhsFilename);
    auto BCPUOperator = std::make_shared<Dune::MatrixAdapter<SpMatrix, Vector, Vector>>(B);

    auto wc = []() -> Vector { throw std::runtime_error("getQuasiImpesWeights is not supported in the benchmarking library."); return Vector(); };


    const size_t N = B.N();

    auto precCPU = PrecFactory::create(*BCPUOperator, configurationCPU.get_child("preconditioner"), wc, 1);
    Dune::InverseOperatorResult resultCPU;

    auto scalarProductCPU = std::make_shared<Dune::SeqScalarProduct<Vector>>();

    auto solverCPU = Solver<Vector>(BCPUOperator,
                                    scalarProductCPU,
                                    precCPU,
                                    configurationCPU.get<double>("tol"),
                                    configurationCPU.get<int>("maxiter"),
                                    configurationCPU.get<int>("verbosity"));

    auto cpuStart = std::chrono::high_resolution_clock::now();
    solverCPU.apply(x, rhs, resultCPU);
    auto cpuEnd = std::chrono::high_resolution_clock::now();

    
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(cpuEnd - cpuStart);
    return std::make_tuple(duration.count(), resultCPU, false);
    
}
template <int dim, template <class> class Solver = Dune::BiCGSTABSolver, class T = double>
std::tuple<unsigned long long, Dune::InverseOperatorResult, bool>
readAndSolveGPU(const auto jsonConfigGPUFilename,
                const auto xFilename,
                const auto matrixFilename,
                const auto rhsFilename)
{
    using M = Opm::MatrixBlock<T, dim, dim>;
    using SpMatrix = Dune::BCRSMatrix<M>;
    using Vector = Dune::BlockVector<Dune::FieldVector<T, dim>>;
    using CuILU0 = Opm::cuistl::CuSeqILU0<SpMatrix, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>;
    using GPUVector = Opm::cuistl::CuVector<T>;
    using Operator = Dune::MatrixAdapter<SpMatrix, Vector, Vector>;
    using PrecFactory = Opm::PreconditionerFactory<Operator, Dune::Amg::SequentialInformation>;

    Opm::PropertyTree configurationGPU(jsonConfigGPUFilename);

    SpMatrix B;
    Vector x, rhs;

    using MatrixDuneIsExpecting = Dune::BCRSMatrix<Dune::FieldMatrix<double, dim, dim>>;

    Dune::loadMatrixMarket(reinterpret_cast<MatrixDuneIsExpecting&>(B), matrixFilename);
    Dune::loadMatrixMarket(x, xFilename);
    Dune::loadMatrixMarket(rhs, rhsFilename);
    auto BCPUOperator = std::make_shared<Dune::MatrixAdapter<SpMatrix, Vector, Vector>>(B);

     auto wc = []() -> Vector { throw std::runtime_error("getQuasiImpesWeights is not supported in the benchmarking library."); return Vector(); };



    const size_t N = B.N();

    Dune::InverseOperatorResult result;

    auto BonGPU = Opm::cuistl::CuSparseMatrix<T>::fromMatrix(B);
    auto BOperator = std::make_shared<
        Dune::MatrixAdapter<Opm::cuistl::CuSparseMatrix<T>, Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>>(
        BonGPU);

    auto precGPUWrapped = PrecFactory::create(*BCPUOperator, configurationGPU.get_child("preconditioner"), wc, 1);

    auto precAsHolder = std::dynamic_pointer_cast<
        Opm::cuistl::PreconditionerHolder<Opm::cuistl::CuVector<T>, Opm::cuistl::CuVector<T>>>(precGPUWrapped);
    if (!precAsHolder) {
        OPM_THROW(std::invalid_argument,
                  "The preconditioner needs to be a CUDA preconditioner wrapped in a "
                  "Opm::cuistl::PreconditionerHolder (eg. CuILU0).");
    }
    auto preconditionerOnGPU = precAsHolder->getUnderlyingPreconditioner();

    auto scalarProduct = std::make_shared<Dune::SeqScalarProduct<Opm::cuistl::CuVector<T>>>();
    auto solver = Solver<Opm::cuistl::CuVector<T>>(BOperator,
                                                   scalarProduct,
                                                   preconditionerOnGPU,
                                                   configurationGPU.get<double>("tol"),
                                                   configurationGPU.get<int>("maxiter"),
                                                   configurationGPU.get<int>("verbosity"));


    auto xOnGPU = GPUVector(x.dim());
    xOnGPU.copyFromHost(x);

    auto rhsOnGPU = GPUVector(rhs.dim());
    rhsOnGPU.copyFromHost(rhs);

    bool gpufailed = false;

    auto gpuStart = std::chrono::high_resolution_clock::now();
    try {
        solver.apply(xOnGPU, rhsOnGPU, result);
    } catch (const std::logic_error& e) {
        gpufailed = true;
    }

    OPM_CUDA_SAFE_CALL(cudaDeviceSynchronize());
    auto gpuEnd = std::chrono::high_resolution_clock::now();

    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(gpuEnd - gpuStart);
    return std::make_tuple(duration.count(), result, gpufailed);
}

boost::property_tree::ptree makeTree(const std::tuple<unsigned long long, Dune::InverseOperatorResult, bool>& resultstuple) {
    boost::property_tree::ptree tree;
    auto [runtime, result, failed] = resultstuple;
    tree.add("runtime", runtime);
    tree.add("failed_by_exception", failed);
    tree.add("iterations", result.iterations);
    tree.add("reduction", result.reduction);
    tree.add("converged", result.converged);
    tree.add("conv_rate", result.conv_rate);
    tree.add("elapsed", result.elapsed);
    tree.add("condition_estimate", result.condition_estimate);

    return tree;
}

template <int dim, class T = double>
void
readAndSolve(const auto jsonConfigCPUFilename,
             const auto jsonConfigGPUFilename,
             const auto xFilename,
             const auto matrixFilename,
             const auto rhsFilename,
             bool runCPU = true,
             bool runGPU = true)
{
    auto tree = boost::property_tree::ptree();

    if (runCPU) {
        Opm::PropertyTree configurationCPU(jsonConfigCPUFilename);
        auto solverNameCPU = configurationCPU.get<std::string>("solver");

        if (solverNameCPU != "bicgstab" && solverNameCPU != "onlypreconditioner") {
            throw std::runtime_error(
                fmt::format("We only support bicgstab and onlypreconditioner, but given {}.", solverNameCPU));
        }

        if (solverNameCPU == "bicgstab") {
            auto runtimeCPU = makeTree(readAndSolveCPU<dim, Dune::BiCGSTABSolver>(
                jsonConfigCPUFilename, xFilename, matrixFilename, rhsFilename));

            tree.add_child("CPU", runtimeCPU);
        } else {
            const auto runtimeCPU = makeTree(readAndSolveCPU<dim, OnlyPreconditionSolver>(
                jsonConfigCPUFilename, xFilename, matrixFilename, rhsFilename));

            tree.add_child("CPU", runtimeCPU);
        }
    }
    if (runGPU) {
        Opm::PropertyTree configurationGPU(jsonConfigGPUFilename);
        auto solverNameGPU = configurationGPU.get<std::string>("solver");

        if (solverNameGPU != "cubicgstab" && solverNameGPU != "onlypreconditioner") {
            throw std::runtime_error(fmt::format(
                "We only support cubicgstab and onlypreconditioner for the GPU, but given {}.", solverNameGPU));
        }

        if (solverNameGPU == "cubicgstab") {
            const auto runtimeGPU = makeTree(readAndSolveGPU<dim, Dune::BiCGSTABSolver>(
                jsonConfigGPUFilename, xFilename, matrixFilename, rhsFilename));
            tree.add_child("GPU", runtimeGPU);

        } else {
            const auto runtimeGPU = makeTree(readAndSolveGPU<dim, OnlyPreconditionSolver>(
                jsonConfigGPUFilename, xFilename, matrixFilename, rhsFilename));
            tree.add_child("GPU", runtimeGPU);
        }
    }


    boost::property_tree::write_json(std::cout, tree, true);
}


int
main(int argc, char** argv)
{

    [[maybe_unused]] const auto& helper = Dune::MPIHelper::instance(argc, argv);

    namespace po = boost::program_options;
    po::options_description desc("Run matrix benchmark.");
    desc.add_options()("help", "Produce this help message.")(
        "matrix-file,m", po::value<std::string>()->required(), "Matrix filename.")(
        "initial-guess-file,x", po::value<std::string>()->required(), "x (initial guess) filename.")(
        "rhs-file,y", po::value<std::string>()->required(), "y (right hand side) filename.")(
        "configfile-cpu,c", po::value<std::string>(), "Configuration file for the linear solver on the CPU (.json)")(
        "configfile-gpu,g", po::value<std::string>(), "Configuration file for the linear solver on the GPU (.json)");

    po::variables_map vm;

    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);

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
    const bool runCPU = vm.count("configfile-cpu") > 0;
    const bool runGPU = vm.count("configfile-gpu") > 0;
    if (!runGPU && !runCPU) {
        std::cerr << "You need to specify at least one of --configfile-cpu or --configfile-gpu\n";
        std::cout << desc << "\n";
        std::exit(EXIT_FAILURE);
    }

    const auto matrixFilename = vm["matrix-file"].as<std::string>();
    const auto xFilename = vm["initial-guess-file"].as<std::string>();
    const auto rhsFilename = vm["rhs-file"].as<std::string>();
    auto jsonConfigFileCPU = std::string("");
    if (runCPU) {
        jsonConfigFileCPU = vm["configfile-cpu"].as<std::string>();
    }

    auto jsonConfigFileGPU = std::string("");
    if (runGPU) {
        jsonConfigFileGPU = vm["configfile-gpu"].as<std::string>();
    }

    size_t dim = 2;
    {
        std::ifstream matrixfile(matrixFilename);
        std::string line;
        const std::string lineToFind = "% ISTL_STRUCT blocked";
        while (std::getline(matrixfile, line)) {
            if (line.substr(0, lineToFind.size()) == lineToFind) {
                dim = std::atoi(line.substr(lineToFind.size() + 2).c_str()); // Yeah, max 9
                break;
            }
        }
    }

    switch (dim) {
    case 1:
        readAndSolve<1>(jsonConfigFileCPU, jsonConfigFileGPU, xFilename, matrixFilename, rhsFilename, runCPU, runGPU);
        break;
    case 2:
        readAndSolve<2>(jsonConfigFileCPU, jsonConfigFileGPU, xFilename, matrixFilename, rhsFilename, runCPU, runGPU);
        break;
    case 3:
        readAndSolve<3>(jsonConfigFileCPU, jsonConfigFileGPU, xFilename, matrixFilename, rhsFilename, runCPU, runGPU);
        break;
    case 4:
        readAndSolve<4>(jsonConfigFileCPU, jsonConfigFileGPU, xFilename, matrixFilename, rhsFilename, runCPU, runGPU);
        break;
    default:
        throw std::runtime_error("Unresolved matrix dimension " + std::to_string(dim));
    }
}
