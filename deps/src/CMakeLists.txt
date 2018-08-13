project(Singular)

cmake_minimum_required(VERSION 3.01)

find_package(JlCxx REQUIRED)

include_directories(${JlCxx_DIR}/../../../include/)
include_directories(${JuliaIncludeDir})
include_directories(${singular_includes})

SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -L${JULIA_LIB_DIR} -Wl,-rpath,${JULIA_LIB_DIR} -L${singular_libdir} -Wl,-rpath,${singular_libdir}" )

add_library(singular SHARED singular.cpp)
target_link_libraries(singular JlCxx::cxxwrap_julia -ljulia "-lSingular -ldl")

install(TARGETS
  singular
ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})