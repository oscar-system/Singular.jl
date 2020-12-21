import Singular_jll, CxxWrap

# TODO: use ARGS to specify custom build dir?
builddir = "build"
installdir = abspath("install")
JlCxx_DIR = joinpath(CxxWrap.prefix_path(), "lib", "cmake", "JlCxx")

rm(builddir; force=true, recursive=true)

run(`cmake
    -DJulia_EXECUTABLE=$(joinpath(Sys.BINDIR, Base.julia_exename()))
    -Dextra_cppflags=-I$(Singular_jll.GMP_jll.artifact_dir)/include
    -Dextra_ldflags=-L$(Singular_jll.GMP_jll.artifact_dir)/lib
    -DSingular_PREFIX=$(Singular_jll.artifact_dir)
    -DCMAKE_INSTALL_PREFIX=$(installdir)
    -DJlCxx_DIR=$(JlCxx_DIR)
    -DCMAKE_BUILD_TYPE=Release
    -S .
    -B $(builddir)
`)

run(`cmake --build $(builddir) --config Release --target install -- -j$(Sys.CPU_THREADS)`)
