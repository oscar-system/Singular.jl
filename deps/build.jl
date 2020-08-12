using BinaryProvider

using Singular_jll

import CxxWrap
import Libdl
import CMake

using GMP_jll
using MPFR_jll

localprefixpath = joinpath(@__DIR__, "usr")

prefixpath = Singular_jll.artifact_dir

gmp_dir = GMP_jll.artifact_dir
mpfr_dir = MPFR_jll.artifact_dir

# add shell scripts that startup another julia for some lib4ti2 programs
# singular and libsingular are supposed to at least look in
#  $prefixpath/lib/singular/MOD

mkpath("$prefixpath/lib/singular/MOD")

write("$prefixpath/lib/singular/MOD/graver", """
#!/bin/sh
julia --startup-file=no -O0 -e 'import lib4ti2_jll
lib4ti2_jll.zsolve() do x
  p=run(ignorestatus(`\$x -G \$ARGS`))
  exit(p.exitcode)
end' -- "\$@"
""")

write("$prefixpath/lib/singular/MOD/hilbert", """
#!/bin/sh
julia --startup-file=no -O0 -e 'import lib4ti2_jll
lib4ti2_jll.zsolve() do x
  p=run(ignorestatus(`\$x -H \$ARGS`))
  exit(p.exitcode)
end' -- "\$@"
""")

write("$prefixpath/lib/singular/MOD/markov", """
#!/bin/sh
julia --startup-file=no -O0 -e 'import lib4ti2_jll
lib4ti2_jll.exe4ti2gmp() do x
  p=run(ignorestatus(`\$x markov \$ARGS`))
  exit(p.exitcode)
end' -- "\$@"
""")

chmod("$prefixpath/lib/singular/MOD/graver", 0o777)
chmod("$prefixpath/lib/singular/MOD/hilbert", 0o777)
chmod("$prefixpath/lib/singular/MOD/markov", 0o777)


push!(Libdl.DL_LOAD_PATH, joinpath(prefixpath, "lib"))

@show libcxxwrap_prefix = CxxWrap.prefix_path()
@show julia_exec = joinpath(Sys.BINDIR, "julia")

extra_cppflags = ""

if Sys.isapple()
   # Work around a bug in Xcode 11.4 that causes SIGABRT, at least until
   # https://github.com/llvm/llvm-project/commit/2464d8135e
   # arrives.
   #
   # We build a jlcxx library that uses std::string which will
   # abort with a failed assertion has_julia_type<T> if we are building
   # with xcode 11.4 but libcxxwrap-julia was built with an older libc++.
   #
   # Read the above LLVM commit message for some details; the effect of merged
   # vs non-merged type_info is that for the former memory addresses are
   # used as hash_code(), for the latter the type_info.name() string is
   # hashed.

   xcodetypeinfo_build_path = joinpath(@__DIR__, "xcodetypeinfo", "build")
   rm(xcodetypeinfo_build_path, recursive = true, force = true)
   mkpath(xcodetypeinfo_build_path)
   cd(xcodetypeinfo_build_path)
   run(`$(CMake.cmake)
        -DJulia_EXECUTABLE=$julia_exec
        -DCMAKE_PREFIX_PATH=$libcxxwrap_prefix
        ..`)
   run(`$(CMake.cmake) --build .`)
   libpath = joinpath(xcodetypeinfo_build_path, "libhello.$(Libdl.dlext)")
   res = run(pipeline(Cmd(`$(Base.julia_cmd()) --project
                           -e "using CxxWrap; @wrapmodule(\"$libpath\", :define_module_hello); @initcxx;"`,
                           ignorestatus = true),
                      stdout=devnull,
                      stderr=devnull))
   if res.termsignal == 6
      println("Applying Xcode type_info.hash_code() workaround")
      global extra_cppflags *= " -DFORCE_XCODE_TYPEINFO_MERGED"
   end
end

cmake_src_path = joinpath(@__DIR__, "src")
cmake_build_path = joinpath(@__DIR__, "build")

rm(cmake_build_path, recursive = true, force = true)
mkpath(cmake_build_path)
cd(cmake_build_path)

@show "$mpfr_dir/lib"
@show "$gmp_dir/lib"

println("Initializing cmake")
run(`$(CMake.cmake)
    -DJulia_EXECUTABLE=$julia_exec
    -DCMAKE_PREFIX_PATH=$libcxxwrap_prefix
    -Dextra_cppflags=$extra_cppflags
    -Dsingular_includes=$prefixpath/include
    -Dsingular_libdir=$prefixpath/lib
    -DCMAKE_INSTALL_LIBDIR=$localprefixpath/lib
    -Dgmp_dir=$gmp_dir
    -Dmpfr_dir=$mpfr_dir
    $cmake_src_path`)

println("Running cmake")
run(`$(CMake.cmake) --build .`)

rm(localprefixpath, recursive = true, force = true)
run(`$(CMake.cmake) --install .`)

include("parselibs.jl")

