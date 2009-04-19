# cmake introduction:
# http://www.linux-magazin.de/heft_abo/ausgaben/2007/02/mal_ausspannen/(offset)/4

cmake_minimum_required(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# project related information
PROJECT(Sassena)
SET(CMAKE_BUILD_TYPE Release)

# some global compiler flags
ADD_DEFINITIONS(-O3 -fPIC -Wall -DBOOST_UBLAS_NDEBUG -DNDEBUG  -funroll-loops -ftree-vectorizer-verbose=3 -pg)
SET(CMAKE_EXE_LINKER_FLAGS -pg -lc_p)


# the elegant way to set up boost
SET(BOOST_ROOT ${EXTERNAL_ROOT})
SET(Boost_ADDITIONAL_VERSIONS "1.37")
# accumulators are included in boost 1.36+
FIND_PACKAGE(Boost 1.36 REQUIRED COMPONENTS regex mpi serialization system filesystem)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}/geometry)

FIND_PACKAGE(MPI REQUIRED)
SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})

# set the libconfig++ library
SET(CONFIG_ROOT ${EXTERNAL_ROOT})
FIND_PACKAGE(Config REQUIRED)
INCLUDE_DIRECTORIES(${CONFIG_INCLUDE_DIR})

# set the global link dependencies and include dirs
SET(LIB_DEPENDENCIES ${CONFIG_LIBRARIES} ${Boost_LIBRARIES})
INCLUDE_DIRECTORIES( ${Sassena_SOURCE_DIR}/include)
MESSAGE(STATUS "sassena_source_dir= " ${Sassena_INCLUDE_DIRS})

ADD_EXECUTABLE(sassena src/sassena.cpp src/sample.cpp src/frame.cpp src/frames.cpp src/atom.cpp src/atoms.cpp src/atomselection.cpp src/analysis.cpp src/coor3d.cpp src/grid.cpp src/settings.cpp src/scatterdata.cpp)

TARGET_LINK_LIBRARIES (sassena ${LIB_DEPENDENCIES})

INSTALL(TARGETS sassena RUNTIME DESTINATION bin)

SET(CPACK_PACKAGE_NAME "Sassena")
SET(CPACK_PACKAGE_VENDOR "Benjamin Lindner (ben@benlabs.net), Franci Merzel (franci.merzel@ki.si)")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Sassena: a method for calculating small-angle X-ray and neutron scattering from explicit-atom models of solvated particles")
SET(CPACK_PACKAGE_VERSION "2.0.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "2")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")

INCLUDE(CPack)

ENABLE_TESTING()
# unit tests:
ADD_EXECUTABLE(unit_coor3d tests/unit_coor3d.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (unit_coor3d ${LIB_DEPENDENCIES})
INSTALL(TARGETS unit_coor3d RUNTIME DESTINATION bin/units)

ADD_EXECUTABLE(unit_boostmpi tests/unit_boostmpi.cpp)
TARGET_LINK_LIBRARIES (unit_boostmpi ${LIB_DEPENDENCIES})
INSTALL(TARGETS unit_boostmpi RUNTIME DESTINATION bin/units)

ADD_EXECUTABLE(unit_exp tests/unit_exp.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (unit_exp ${LIB_DEPENDENCIES})
INSTALL(TARGETS unit_exp RUNTIME DESTINATION bin/units)

ADD_EXECUTABLE(unit_libconfig tests/unit_libconfig.cpp)
TARGET_LINK_LIBRARIES (unit_libconfig ${LIB_DEPENDENCIES})
INSTALL(TARGETS unit_libconfig RUNTIME DESTINATION bin/units)

ADD_EXECUTABLE(unit_triangulation tests/unit_triangulation.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (unit_triangulation ${LIB_DEPENDENCIES})
INSTALL(TARGETS unit_triangulation RUNTIME DESTINATION bin/units)

ADD_EXECUTABLE(unit_randomnumbers tests/unit_randomnumbers.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (unit_randomnumbers ${LIB_DEPENDENCIES})
INSTALL(TARGETS unit_randomnumbers RUNTIME DESTINATION bin/units)

# additional tools:
ADD_EXECUTABLE(normalize tools/normalize.cpp)
TARGET_LINK_LIBRARIES (normalize )
INSTALL(TARGETS normalize RUNTIME DESTINATION bin/tools)

ADD_EXECUTABLE(grid tools/grid.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (grid )
INSTALL(TARGETS grid RUNTIME DESTINATION bin/tools)


# additional tools:
ADD_EXECUTABLE(map tools/map.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (map )
INSTALL(TARGETS map RUNTIME DESTINATION bin/tools)

ADD_TEST(unit_coor3d "Coordinate 3D type system")
ADD_TEST(unit_exp "Boost Special functions")
ADD_TEST(unit_boostmpi  "Boost MPI implementation")
ADD_TEST(unit_libconfig  "Libconfig implementation")
