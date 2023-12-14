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
find_package(amgcl)

if (amgcl_FOUND)
    SET(HAVE_AMGCL ON)
else()
    SET(HAVE_AMGCL OFF)
endif()

if (OPM_REQUIRES_DAMARIS)
  find_package(XercesC REQUIRED)
  find_package(Damaris REQUIRED)
endif()

if (OPM_REQUIRES_VEXCL)
  include(VexCLSetup)
endif()

find_package(Boost REQUIRED COMPONENTS program_options thread log log_setup python numpy unit_test_framework)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/config.hpp ESCAPE_QUOTES)

add_executable(linsolverlab src/main.cpp
    src/read_binary.hpp)

target_compile_features(linsolverlab PRIVATE cxx_std_20)
target_include_directories(linsolverlab PRIVATE 
  # Does not seem like these are incluced in any target? Doing it the old fashioned way...
  ${opm-simulators_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}/include
)

target_link_libraries(linsolverlab PRIVATE opmsimulators Boost::program_options CUDA::cudart fmt::fmt)

if (HAVE_AMGCL)
    target_link_libraries(linsolverlab PRIVATE opmsimulators amgcl::amgcl)


    add_executable(test_read_binary tests/test_read_binary.cpp)
    target_compile_features(test_read_binary PRIVATE cxx_std_20)

    target_include_directories(test_read_binary PRIVATE
      # Does not seem like these are incluced in any target? Doing it the old fashioned way...
      ${opm-simulators_INCLUDE_DIRS}
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_CURRENT_BINARY_DIR}/include
    )

    target_link_libraries(test_read_binary PRIVATE opmsimulators Boost::program_options Boost::unit_test_framework amgcl::amgcl CUDA::cudart fmt::fmt)

endif()

if(Petsc_FOUND)
  add_executable(linearsolverlab_petsc src/petsc.cpp)
  target_include_directories(linearsolverlab_petsc PRIVATE ${PETSC_INCLUDE_DIRS} ${opm-simulators_INCLUDE_DIRS})
  target_link_libraries(linearsolverlab_petsc PRIVATE ${PETSC_LIBRARIES} fmt::fmt Boost::program_options ${opm-simulators_LIBRARIES})
endif()

