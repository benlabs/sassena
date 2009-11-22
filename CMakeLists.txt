# cmake introduction:
# http://www.linux-magazin.de/heft_abo/ausgaben/2007/02/mal_ausspannen/(offset)/4

cmake_minimum_required(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules")

# project related information
PROJECT(Sassena)
SET(CMAKE_BUILD_TYPE Release)

## some global compiler flags
SET(WALLFLAGS ${WALLFLAGS} "-ftree-vectorize  ")
##SET(WALLFLAGS ${WALLFLAGS} "-Waddress                        "        )
##SET(WALLFLAGS ${WALLFLAGS} "-Warray-bounds                   "        )
##SET(WALLFLAGS ${WALLFLAGS} "-Wc++0x-compat                   "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wchar-subscripts                "        )
##SET(WALLFLAGS ${WALLFLAGS} "-Wenum-compare              		"  )
##SET(WALLFLAGS ${WALLFLAGS} "-Wimplicit-int                   "        )
##SET(WALLFLAGS ${WALLFLAGS} "-Wimplicit-function-declaration  "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wcomment                        "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wformat                         "        )
##SET(WALLFLAGS ${WALLFLAGS} "-Wmain                      		"  )
#SET(WALLFLAGS ${WALLFLAGS} "-Wmissing-braces                 "        )
##SET(WALLFLAGS ${WALLFLAGS} "-Wnonnull                        "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wparentheses                    "        )
##SET(WALLFLAGS ${WALLFLAGS} "-Wpointer-sign                   "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wreorder                        "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wreturn-type                    "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wsequence-point                 "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wsign-compare                   "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wstrict-aliasing                "        )
##SET(WALLFLAGS ${WALLFLAGS} "-Wstrict-overflow              "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wswitch                         "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wtrigraphs                      "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wuninitialized                  "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wunknown-pragmas                "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wunused-function                "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wunused-label                   "        )
#SET(WALLFLAGS ${WALLFLAGS} "-Wunused-value                   "        )
##SET(WALLFLAGS ${WALLFLAGS} -Wunused-variable                   ) 
##SET(WALLFLAGS ${WALLFLAGS} "-Wvolatile-register-var"              )

#ADD_DEFINITIONS( -O3 -fPIC -DBOOST_UBLAS_NDEBUG -DNDEBUG  -funroll-loops -ftree-vectorizer-verbose=3 -pg )
#ADD_DEFINITIONS( -O3 -fPIC )
ADD_DEFINITIONS( -fPIC )
#SET(CMAKE_EXE_LINKER_FLAGS -pg -lc_p)



# the elegant way to set up boost
SET(BOOST_ROOT ${EXTERNAL_ROOT})
SET(Boost_ADDITIONAL_VERSIONS "1.37")
# accumulators are included in boost 1.36+
FIND_PACKAGE(Boost 1.36 REQUIRED COMPONENTS regex mpi serialization system filesystem)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS}/geometry)


FIND_PACKAGE(FFTW3 REQUIRED)
INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_DIR})

FIND_PACKAGE(MPI REQUIRED)
SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})

# set the libconfig++ library
SET(CONFIG_ROOT ${EXTERNAL_ROOT})
FIND_PACKAGE(Config REQUIRED)
INCLUDE_DIRECTORIES(${Sassena_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CONFIG_INCLUDE_DIR})

# libxml support
FIND_PACKAGE(LibXml2 REQUIRED)
INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})

# set the global link dependencies and include dirs
SET(LIB_DEPENDENCIES ${CONFIG_LIBRARIES} ${Boost_LIBRARIES} ${LIBXML2_LIBRARIES} ${FFTW3_LIBRARIES} ) 
INCLUDE_DIRECTORIES( ${Sassena_SOURCE_DIR}/include/xdrfile)


ADD_LIBRARY(xdrfile-static STATIC src/xdrfile/xdrfile.c src/xdrfile/xdrfile_xtc.c src/xdrfile/xdrfile_trr.c )
SET_TARGET_PROPERTIES(xdrfile-static PROPERTIES OUTPUT_NAME "xdrfile")
SET(LIB_DEPENDENCIES ${LIB_DEPENDENCIES} xdrfile-static)

MESSAGE(STATUS "sassena_lib_depend= " ${LIB_DEPENDENCIES})

ADD_LIBRARY(sass_control SHARED
	src/control/log.cpp
	src/control/parameters.cpp
	src/control/database.cpp
	src/xml_interface.cpp
	src/coor3d.cpp
)

TARGET_LINK_LIBRARIES (sass_control ${LIB_DEPENDENCIES} )
SET(LIB_DEPENDENCIES ${LIB_DEPENDENCIES} sass_control)


ADD_LIBRARY(sass_sample SHARED
	src/sample/sample.cpp
	src/sample/atom.cpp
	src/sample/atoms.cpp
	src/sample/atomselection.cpp
	src/sample/coordinate_set.cpp
	src/sample/coordinate_sets.cpp
	src/sample/frame.cpp
	src/sample/frames.cpp
	src/sample/motion_walker.cpp	
	
	src/measures/center_of_mass.cpp
	src/decompose.cpp
)

TARGET_LINK_LIBRARIES (sass_sample ${LIB_DEPENDENCIES} )
SET(LIB_DEPENDENCIES ${LIB_DEPENDENCIES} sass_sample)


ADD_EXECUTABLE(sassena 
    src/main/sassena.cpp

	src/analysis.cpp
	src/common.cpp
	src/coor3d.cpp 
	src/decompose.cpp
	src/decomposition_plan.cpp
	src/density_grid.cpp
	src/fftw/fftw++.cc
	src/grid.cpp  
	src/measures/center_of_mass.cpp
	src/performance_analyzer.cpp
	src/progress_reporter.cpp
	src/scatter_devices/all_exact.cpp 	
	src/scatter_devices/all_multipole.cpp 	
	src/scatter_devices/all_vectors.cpp 
	src/scatter_devices/self_vectors.cpp 
	src/scatter_factors.cpp
	src/scatter_spectrum.cpp 
	src/scatter_spectrum_writer.cpp 
	src/timer.cpp 
	src/xml_interface.cpp
	src/vector_unfold.cpp
	)
	
	
ADD_EXECUTABLE(s_coordump 
    src/main/s_coordump.cpp

	src/common.cpp
	src/coor3d.cpp 
	src/decompose.cpp
	src/measures/center_of_mass.cpp
	src/xml_interface.cpp
	src/vector_unfold.cpp
	)


	
ADD_EXECUTABLE(s_savesample 
    src/main/s_savesample.cpp

	src/common.cpp
	src/coor3d.cpp 
	src/decompose.cpp
	src/measures/center_of_mass.cpp
	src/xml_interface.cpp
	src/vector_unfold.cpp
	)
	
TARGET_LINK_LIBRARIES (sassena ${LIB_DEPENDENCIES} )
TARGET_LINK_LIBRARIES (s_coordump ${LIB_DEPENDENCIES} )
TARGET_LINK_LIBRARIES (s_savesample ${LIB_DEPENDENCIES} )

INSTALL(TARGETS sassena RUNTIME DESTINATION bin)
INSTALL(TARGETS s_coordump RUNTIME DESTINATION bin)
INSTALL(TARGETS s_savesample RUNTIME DESTINATION bin)

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
#ADD_EXECUTABLE(unit_coor3d tests/unit_coor3d.cpp src/coor3d.cpp)
#TARGET_LINK_LIBRARIES (unit_coor3d ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_coor3d RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_boostmpi tests/unit_boostmpi.cpp)
#TARGET_LINK_LIBRARIES (unit_boostmpi ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_boostmpi RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_exp tests/unit_exp.cpp src/coor3d.cpp)
#TARGET_LINK_LIBRARIES (unit_exp ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_exp RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_libconfig tests/unit_libconfig.cpp)
#TARGET_LINK_LIBRARIES (unit_libconfig ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_libconfig RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_triangulation tests/unit_triangulation.cpp src/coor3d.cpp)
#TARGET_LINK_LIBRARIES (unit_triangulation ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_triangulation RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_randomnumbers tests/unit_randomnumbers.cpp src/coor3d.cpp)
#TARGET_LINK_LIBRARIES (unit_randomnumbers ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_randomnumbers RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_ranklogic tests/unit_ranklogic.cpp src/tasks.cpp src/coor3d.cpp)
#TARGET_LINK_LIBRARIES (unit_ranklogic ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_ranklogic RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_serializestatic tests/unit_serializestatic.cpp )
#TARGET_LINK_LIBRARIES (unit_serializestatic ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_serializestatic RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_mpicollectives tests/unit_mpicollectives.cpp src/timer.cpp src/log.cpp)
#TARGET_LINK_LIBRARIES (unit_mpicollectives ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_mpicollectives RUNTIME DESTINATION bin/units)
#
#ADD_EXECUTABLE(unit_sasstest 
#	tests/unit_sasstest.cpp 
#	src/sample.cpp src/parameters.cpp src/log.cpp src/frame.cpp src/frames.cpp src/atom.cpp src/atoms.cpp src/atomselection.cpp src/coor3d.cpp src/settings.cpp src/tasks.cpp src/decompose.cpp src/scatterdevices.cpp src/coordinateset.cpp src/scatterspectrum.cpp src/analysis.cpp src/grid.cpp)
#TARGET_LINK_LIBRARIES (unit_sasstest ${LIB_DEPENDENCIES})
#INSTALL(TARGETS unit_sasstest RUNTIME DESTINATION bin/units)


# additional tools:
ADD_EXECUTABLE(normalize tools/normalize.cpp)
TARGET_LINK_LIBRARIES (normalize )
INSTALL(TARGETS normalize RUNTIME DESTINATION bin/tools)

ADD_EXECUTABLE(grid tools/grid.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (grid )
INSTALL(TARGETS grid RUNTIME DESTINATION bin/tools)

ADD_EXECUTABLE(grid_frames tools/grid_frames.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (grid_frames )
INSTALL(TARGETS grid_frames RUNTIME DESTINATION bin/tools)


# additional tools:
ADD_EXECUTABLE(map tools/map.cpp src/coor3d.cpp)
TARGET_LINK_LIBRARIES (map )
INSTALL(TARGETS map RUNTIME DESTINATION bin/tools)

ADD_TEST(unit_coor3d "Coordinate 3D type system")
ADD_TEST(unit_exp "Boost Special functions")
ADD_TEST(unit_boostmpi  "Boost MPI implementation")
ADD_TEST(unit_libconfig  "Libconfig implementation")
