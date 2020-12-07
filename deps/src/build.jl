#cmake -DJulia_EXECUTABLE=/home/bla/Downloads/julia-1.5.3/bin/julia -DSingular_PREFIX=/home/bla/.julia/dev/Singular/deps/usr/ -DCMAKE_INSTALL_PREFIX=/tmp -DJlCxx_DIR=/home/bla/.julia/artifacts/860a8b2216bd059600ed7c44cdaa3bb81b23ff1c/lib/cmake/JlCxx -DCMAKE_BUILD_TYPE=Release -S libsingular-julia -B /tmp/build

import Singular, CxxWrap

Julia_EXECUTABLE = joinpath(Sys.BINDIR, Base.julia_exename())
Singular_PREFIX = joinpath(dirname(pathof(Singular)), "..", "deps", "usr")
JlCxx_DIR = joinpath(CxxWrap.CxxWrapCore.libcxxwrap_julia_jll.artifact_dir, "lib", "cmake", "JlCxx")

cm = `cmake -DJulia_EXECUTABLE=$(Julia_EXECUTABLE) -DSingular_PREFIX=$(Singular_PREFIX) -DCMAKE_INSTALL_PREFIX=/tmp -DJlCxx_DIR=$(JlCxx_DIR) -DCMAKE_BUILD_TYPE=Release -S libsingular-julia -B /tmp/build`

run(cm)

run(`cmake --build /tmp/build --config Release --target install`)
