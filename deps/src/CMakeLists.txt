cmake_minimum_required(VERSION 3.01)

# I don't know why I need the following, but it works
if(POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)
endif()

project(libsingular_julia)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(JlCxx REQUIRED)

include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 14)

set(JLSINGULAR_TARGET singular_julia)

# avoid gcc 9 internal compiler error,
# see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90998
if(CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.0
    AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.3)
  SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wno-conversion" )
endif()

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -std=c++14" )

# to avoid lib64 dirs as binarybuilder uses lib everywhere
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_INSTALL_LIBDIR "lib")

get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

include_directories(${CMAKE_INSTALL_PREFIX})
include_directories(${Singular_PREFIX}/include/singular)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g ${extra_cppflags} -I${CMAKE_INSTALL_PREFIX}/include -I${Singular_PREFIX}/include/ " )
SET( CMAKE_SHARED_LINKER_FLAGS  "${CMAKE_SHARED_LINKER_FLAGS} -v -L${CMAKE_INSTALL_PREFIX}/lib -Wl,-rpath,${CMAKE_INSTALL_PREFIX}/lib -L${Singular_PREFIX}/lib -Wl,-rpath,${Singular_PREFIX}/lib" )


add_library(${JLSINGULAR_TARGET} SHARED singular.cpp rings.cpp coeffs.cpp ideals.cpp matrices.cpp caller.cpp coeff_rings.cpp threading.cpp)

target_link_libraries(${JLSINGULAR_TARGET} JlCxx::cxxwrap_julia -lSingular -lpolys -lsingular_resources -lfactory -lomalloc -ldl -lgmp)

install(TARGETS
  ${JLSINGULAR_TARGET}
ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
