project(Singular)

cmake_minimum_required(VERSION 3.01)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(JlCxx REQUIRED)

get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

include_directories(${singular_includes})
include_directories(${singular_includes}/singular)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ${extra_cppflags} -I${gmp_dir}/include -I${mpfr_dir}/include" )
SET( CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -v -L${singular_libdir} -L${gmp_dir}/lib -L${mpfr_dir}/lib -Wl,-rpath,${singular_libdir} -Wl,-rpath,${gmp_dir}/lib -Wl,-rpath,${mpfr_dir}/lib" )


add_library(singularwrap SHARED singular.cpp rings.cpp coeffs.cpp ideals.cpp matrices.cpp caller.cpp coeff_rings.cpp threading.cpp)

target_link_libraries(singularwrap JlCxx::cxxwrap_julia -lSingular -lpolys -lsingular_resources -lfactory -lomalloc -ldl -lgmp)

install(TARGETS
  singularwrap
ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
