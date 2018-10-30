import Libdl

const oldwdir = pwd()
const pkgdir = dirname(dirname(Base.find_package("Singular")))
const nemodir = dirname(dirname(Base.find_package("Nemo")))

const debug_build = true # N.B: debug builds are up to 50 times slower at runtime!

wdir = "$pkgdir/deps"
vdir = "$pkgdir/local"
nemovdir = "$nemodir/local"

LDFLAGS = "-rpath $vdir/lib -R$vdir/lib -R$nemovdir/lib -R\$\$ORIGIN/../share/julia/site/v$(VERSION.major).$(VERSION.minor)/Singular/local/lib"

cd(wdir)

# INSTALL NTL

const ntl="ntl-9.3.0"

try
  run(`wget -q -nc -c -O "$wdir/$ntl.tar.gz" "http://www.shoup.net/ntl/$ntl.tar.gz"`)
catch
end

const tmp = mktempdir(wdir)

# http://www.shoup.net/ntl/WinNTL-9_3_0.zip
cd(tmp)
run(`tar -C "$tmp" -xkvf "$wdir/$ntl.tar.gz"`)
cd(joinpath(tmp, ntl, "src"))
withenv("CPP_FLAGS"=>"-I$vdir/include", "LD_LIBRARY_PATH"=>"$vdir/lib:$nemovdir/lib", "LDFLAGS"=>LDFLAGS) do
   run(`./configure PREFIX=$vdir DEF_PREFIX=$nemovdir SHARED=on NTL_THREADS=off NTL_EXCEPTIONS=off NTL_GMP_LIP=on CXXFLAGS="-I$nemovdir/include"`)
   run(`make -j4`)
   run(`make install`)
end
cd(wdir)
rm(tmp; recursive=true)

# Install cddlib

# Currently cddlib appears to have no way to specify what GMP is used
# thus --enable-gfanlib is not possible below, since it relies on cddlib
# being installed

# const cddlib="cddlib-094h"

# try
#   run(`wget -q -nc -c -O $wdir/$cddlib.tar.gz ftp://ftp.math.ethz.ch/users/fukudak/cdd/$cddlib.tar.gz`)
# catch
# end

# const tmp2 = mktempdir(wdir)

# cd(tmp2)
# run(`tar -C $tmp2 -xkvf $wdir/$cddlib.tar.gz`)
# cd(joinpath(tmp2, cddlib))
# withenv("CPP_FLAGS"=>"-I$vdir/include", "LD_LIBRARY_PATH"=>"$vdir/lib:$nemovdir/lib", "LDFLAGS"=>LDFLAGS) do
#    run(`./configure --prefix=$vdir --with-gmp=$nemovdir`)
#    run(`make -j4`)
#    run(`make install`)
# end
# cd(wdir)
# rm(tmp2; recursive=true)

# Install Singular

const srcs = joinpath(wdir, "Singular")

# get/update sources
try
  run(`git clone -b spielwiese https://github.com/Singular/Sources.git $srcs`)
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
withenv("CPP_FLAGS"=>"-I$vdir/include", "LD_LIBRARY_PATH"=>"$vdir/lib:$nemodir/lib") do
   if !debug_build
      run(`$srcs/configure --prefix=$vdir --disable-static --disable-p-procs-static --disable-gfanlib --enable-p-procs-dynamic --enable-shared --with-gmp=$nemovdir --with-flint=$nemovdir --with-ntl=$vdir --without-python --with-readline=no`)
   else
      run(`$srcs/configure --prefix=$vdir --disable-static --disable-p-procs-static --disable-gfanlib --enable-p-procs-dynamic --enable-shared --with-gmp=$nemovdir --with-flint=$nemovdir --with-ntl=$vdir --without-python --with-readline=no --with-debug --enable-debug --disable-optimizationflags`)
   end
   withenv("LDFLAGS"=>LDFLAGS) do
      run(`make -j4`)
      run(`make install`)
   end
end

cd(wdir)

push!(Libdl.DL_LOAD_PATH, "$pkgdir/local/lib")

cd(oldwdir)

jlcxx_cmake_dir = abspath(joinpath(dirname(Base.find_package("CxxWrap")), "..", "deps", "usr", "lib", "cmake", "JlCxx"))

julia_include = abspath(joinpath(Sys.BINDIR, "..", "include", "julia"))
julia_lib = abspath(joinpath(Sys.BINDIR, "..", "lib"))
julia_exec = joinpath(Sys.BINDIR, "julia")

cmake_build_path = abspath(joinpath(dirname(Base.find_package("Singular")), "..", "deps", "src"))

cd(cmake_build_path)

run(`cmake -DJulia_EXECUTABLE=$julia_exec -DJlCxx_DIR=$jlcxx_cmake_dir -DJuliaIncludeDir=$julia_include -DJULIA_LIB_DIR=$julia_lib -Dnemo_includes=$nemovdir/include -Dsingular_includes=$vdir/include -Dsingular_libdir=$vdir/lib -DCMAKE_INSTALL_LIBDIR=$vdir/lib .`)

run(`make`)
run(`make install`)

