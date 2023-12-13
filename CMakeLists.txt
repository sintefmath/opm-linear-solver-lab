cmake_minimum_required(VERSION 3.27)
# We need this for Damaris
cmake_policy(SET CMP0148 OLD)
project(
  opm-linear-solver-lab
  VERSION 1.0
  LANGUAGES C CXX CUDA)
option(OPM_REQUIRES_DAMARIS OFF)
option(OPM_REQUIRES_VEXCL OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

find_package(CUDAToolkit REQUIRED)
find_package(TBB REQUIRED)
find_package(dune-common REQUIRED)
find_package(dune-grid REQUIRED)
find_package(dune-geometry REQUIRED)
find_package(dune-istl REQUIRED)
find_package(opm-simulators REQUIRED)
find_package(fmt REQUIRED)
find_package(Petsc)

if (${OPM_REQUIRES_DAMARIS})
  find_package(XercesC REQUIRED)
  find_package(Damaris REQUIRED)
endif()

if (${OPM_REQUIRES_VEXCL})
  include(VexCLSetup)
endif()

find_package(Boost REQUIRED COMPONENTS program_options thread log log_setup python numpy)

add_executable(linsolverlab src/main.cpp)
target_compile_features(linsolverlab PRIVATE cxx_std_20)
target_include_directories(linsolverlab PRIVATE 
  # Does not seem like these are incluced in any target? Doing it the old fashioned way...
  ${opm-simulators_INCLUDE_DIRS}
)


target_link_libraries(linsolverlab PRIVATE opmsimulators Boost::program_options CUDA::cudart fmt::fmt)

if(Petsc_FOUND)
  add_executable(linearsolverlab_petsc src/petsc.cpp)
  target_include_directories(linearsolverlab_petsc PRIVATE ${PETSC_INCLUDE_DIRS} ${opm-simulators_INCLUDE_DIRS})
  target_link_libraries(linearsolverlab_petsc PRIVATE ${PETSC_LIBRARIES} fmt::fmt Boost::program_options ${opm-simulators_LIBRARIES})
endif()
