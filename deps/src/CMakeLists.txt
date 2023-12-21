cmake_minimum_required(VERSION 3.5)

project(libsingular_julia)

set(JL_TARGET singular_julia)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
message(STATUS "Found JlCxx at ${JlCxx_location}")

include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 14)

# avoid gcc 9 internal compiler error,
# see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90998
if(CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.0
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.3)
  SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-conversion" )
endif()

# to avoid lib64 dirs as binarybuilder uses lib everywhere
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_INSTALL_LIBDIR "lib")

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

include_directories(
    ${Singular_PREFIX}/include
    ${Singular_PREFIX}/include/singular
    ${gmp_prefix}/include
    )

link_directories(
    ${Singular_PREFIX}/lib
    ${gmp_prefix}/lib
    )

#SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ${extra_cppflags}" )
#SET( CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -v ${extra_ldflags} )


add_library(${JL_TARGET} SHARED
    singular.cpp
    rings.cpp
    coeffs.cpp
    ideals.cpp
    matrices.cpp
    caller.cpp
    coeff_rings.cpp
    threading.cpp
    )

target_link_libraries(${JL_TARGET}
    JlCxx::cxxwrap_julia
    -lSingular
    -lpolys
    -lsingular_resources
    -lfactory
    -lomalloc
    -ldl
    -lgmp
    )

install(TARGETS
  ${JL_TARGET}
ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
