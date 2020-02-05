using BinaryProvider

import CxxWrap
import Libdl
import Nemo

# Parse some basic command-line arguments
const verbose = "--verbose" in ARGS

issource_build = "SINGULAR_SOURCE_BUILD" in keys(ENV) && ENV["SINGULAR_SOURCE_BUILD"] == "1"

const prefixpath = joinpath(@__DIR__, "usr")
const wdir = joinpath(@__DIR__)
const oldwdir = pwd()

const nemodir = realpath(joinpath(dirname(pathof(Nemo)), ".."))
nemovdir = "$nemodir/deps/usr"

if !issource_build
   # Dependencies that must be installed before this package can be built
   dependencies = [
      # This has to be in sync with the corresponding commit in the source build below (for flint, singular)
      "https://github.com/JuliaPackaging/Yggdrasil/releases/download/GMP-v6.1.2-1/build_GMP.v6.1.2.jl",
      "https://github.com/JuliaPackaging/Yggdrasil/releases/download/MPFR-v4.0.2-1/build_MPFR.v4.0.2.jl",
      "https://github.com/thofma/Flint2Builder/releases/download/ba0cee/build_libflint.v0.0.0-ba0ceed35136a2a43441ab9a9b2e7764e38548ea.jl",
      "https://github.com/thofma/NTLBuilder2/releases/download/v10.5.0-1/build_libntl.v10.5.0.jl",
      "https://github.com/wbhart/SingularBuilder/releases/download/0.0.2/build_libsingular.v0.0.1-2.jl",
   ]

   const prefix = Prefix(get([a for a in ARGS if a != "--verbose"], 1, joinpath(@__DIR__, "usr")))

   products = []

   for url in dependencies
      build_file = joinpath(@__DIR__, basename(url))
      if !isfile(build_file)
         download(url, build_file)
      end
   end

   # Execute the build scripts for the dependencies in an isolated module to
   # avoid overwriting any variables/constants here
   for url in dependencies
      build_file = joinpath(@__DIR__, basename(url))
      m = @eval module $(gensym()); include($build_file); end
      append!(products, m.products)
   end

   filenames = ["libgmp.la", "libgmpxx.la", "libmpfr.la"]

   for filename in filenames
      fpath = joinpath(prefixpath, "lib", filename)
      txt = read(fpath, String)
      open(fpath, "w") do f
         write(f, replace(txt, "/workspace/destdir" => prefixpath))
      end
   end

else # source build

   ("NEMO_SOURCE_BUILD" in keys(ENV) && ENV["NEMO_SOURCE_BUILD"] == "1") || error("Source build of Nemo required")

   const pkgdir = realpath(joinpath(@__DIR__, ".."))

   const debug_build = false # N.B: debug builds are up to 50 times slower at runtime!

   @show NTL_VERSION = "10.5.0"
   @show SINGULAR_VERSION = "d2a9865cca207801d3b0c05a50dfc5a02724120e"
   @show CDDLIB_VERSION = "094h"

   println("Removing old binaries ...")

   rm(prefixpath, force = true, recursive = true)
   mkdir(prefixpath)
   mkdir(joinpath(prefixpath, "lib"))

   if Sys.isapple() && !("CC" in keys(ENV))
      ENV["CC"] = "clang"
      ENV["CXX"] = "clang++"
   end

   LDFLAGS = "-Wl,-rpath,$prefixpath/lib -Wl,-rpath,\$\$ORIGIN/../share/julia/site/v$(VERSION.major).$(VERSION.minor)/Singular/local/lib"
   cd(wdir)

   cores = Sys.CPU_THREADS
   println("Detected $cores CPU threads.")

   # INSTALL NTL

   const NTL_FILE=joinpath(wdir, "ntl-$NTL_VERSION"*".tar.gz")

   try
      download("https://www.shoup.net/ntl/$NTL_FILE", "$NTL_FILE")
   catch
   end

   const tmp = mktempdir(wdir)

   # http://www.shoup.net/ntl/
   cd(tmp)
   run(`tar -C "$tmp" -zxkvf $NTL_FILE`)
   cd(joinpath(tmp, "ntl-$NTL_VERSION", "src"))
   withenv("CPP_FLAGS"=>"-I$prefixpath/include", "LD_LIBRARY_PATH"=>"$prefixpath/lib:$nemovdir/lib", "LDFLAGS"=>LDFLAGS) do
      run(`./configure PREFIX=$prefixpath DEF_PREFIX=$nemovdir SHARED=on NTL_THREADS=off NTL_EXCEPTIONS=off NTL_GMP_LIP=on CXXFLAGS="-I$nemovdir/include"`)
      run(`make -j$cores`)
      run(`make install`)
   end
   cd(wdir)
   rm(tmp; recursive=true)

   # Install cddlib

   # Currently cddlib appears to have no way to specify what GMP is used
   # thus --enable-gfanlib is not possible below, since it relies on cddlib
   # being installed

   # const CDDLIB_FILE="cddlib-$CDDLIB_VERSION"*".tar.gz"

   # try
   #   run(`wget -q -nc -c -O $wdir/$CDDLIB_FILE ftp://ftp.math.ethz.ch/users/fukudak/cdd/$CDDLIB_FILE`)
   # catch
   # end

   # const tmp2 = mktempdir(wdir)

   # cd(tmp2)
   # run(`tar -C $tmp2 -xkvf $wdir/$CDDLIB_FILE`)
   # cd(joinpath(tmp2, cddlib))
   # withenv("CPP_FLAGS"=>"-I$prefixpath/include", "LD_LIBRARY_PATH"=>"$prefixpath/lib:$nemovdir/lib", "LDFLAGS"=>LDFLAGS) do
   #    run(`./configure --prefix=$prefixpath --with-gmp=$nemovdir`)
   #    run(`make -j$cores`)
   #    run(`make install`)
   # end
   # cd(wdir)
   # rm(tmp2; recursive=true)

   # Install Singular

   const srcs = joinpath(wdir, "Singular")

   # get/update sources
   try
      run(`git clone -b spielwiese https://github.com/Singular/Sources.git $srcs`)
      run(`git checkout $SINGULAR_VERSION`)
   catch
      cd(srcs)
      try
         run(`git pull --rebase`)
      catch
      end
      cd(wdir)
   end  

   cd(srcs)
   run(`./autogen.sh`)
   cd(wdir)

   # out of source-tree building:
   try 
      mkdir(joinpath(wdir, "Singular_build"))
   catch
   end

   cd(joinpath(wdir, "Singular_build"))
   withenv("CPP_FLAGS"=>"-I$prefixpath/include", "LD_LIBRARY_PATH"=>"$prefixpath/lib:$nemodir/lib") do
      cmd = split(
        """
        $srcs/configure
        --with-libparse
        --prefix=$prefixpath
        --libdir=$prefixpath/lib
        --disable-static
        --enable-p-procs-static
        --disable-p-procs-dynamic
        --disable-gfanlib
        --enable-shared
        --with-gmp=$nemovdir
        --with-flint=$nemovdir
        --with-ntl=$prefixpath
        --without-python
        --with-readline=no
        """, "\n", keepempty = false)
      if debug_build
         append!(cmd, [
           "--with-debug",
           "--enable-debug",
           "--disable-optimizationflags",
         ])
      end
      run(Cmd(string.(cmd)))
      withenv("LDFLAGS"=>LDFLAGS) do
         run(`make -j$cores`)
         run(`make install`)
      end
   end

   print("Done building Singular")

   cd(wdir)
end

push!(Libdl.DL_LOAD_PATH, joinpath(prefixpath, "lib"))

cd(oldwdir)

println("here")

# jlcxx_cmake_dir = realpath(joinpath(dirname(pathof(CxxWrap)), "..", "deps", "usr", "lib", "cmake", "JlCxx"))

jlcxx_cmake_dir = ENV["JLCXX_DIR"]

println("jlcxx_cmake_dir = $jlcxx_cmake_dir")

julia_include = realpath(joinpath(Sys.BINDIR, "..", "include", "julia"))
julia_lib = realpath(joinpath(Sys.BINDIR, "..", "lib"))
julia_exec = joinpath(Sys.BINDIR, "julia")

cmake_build_path = joinpath(@__DIR__, "src")

cd(cmake_build_path)

print("Initializing cmake")

run(`cmake -DJulia_EXECUTABLE=$julia_exec -DJlCxx_DIR=$jlcxx_cmake_dir -DJuliaIncludeDir=$julia_include -DJULIA_LIB_DIR=$julia_lib -Dnemo_includes=$nemovdir/include -Dsingular_includes=$prefixpath/include -Dsingular_libdir=$prefixpath/lib -DCMAKE_INSTALL_LIBDIR=$prefixpath/lib .`)

print("Running cmake")

run(`make VERBOSE=1`)
run(`make install`)

include("parselibs.jl")

