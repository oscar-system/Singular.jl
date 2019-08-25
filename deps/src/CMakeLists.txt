project(Singular)

cmake_minimum_required(VERSION 3.01)

find_package(JlCxx REQUIRED)

include_directories(${JlCxx_DIR}/../../../include/)
include_directories(${JuliaIncludeDir})
include_directories(${nemo_includes})
include_directories(${singular_includes})
include_directories(${singular_includes}/singular)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14" )
SET( CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -L${JULIA_LIB_DIR} -Wl,-rpath,${JULIA_LIB_DIR} -L${singular_libdir} -Wl,-rpath,${singular_libdir}" )

add_library(singularwrap SHARED singular.cpp rings.cpp coeffs.cpp ideals.cpp matrices.cpp coeff_rings.cpp)
target_link_libraries(singularwrap JlCxx::cxxwrap_julia -ljulia "-lSingular -lpolys -lsingular_resources -lfactory -lomalloc -ldl")

install(TARGETS
  singularwrap
ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
