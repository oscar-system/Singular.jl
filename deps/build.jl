const oldwdir = pwd()
const pkgdir = Pkg.dir("Singular") 
const nemodir = Pkg.dir("Nemo")

const debug_build = false # N.B: debug builds are up to 50 times slower at runtime!

wdir = "$pkgdir/deps"
vdir = "$pkgdir/local"
nemovdir = "$nemodir/local"

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
run(`./configure DEF_PREFIX="$vdir" SHARED=on NTL_THREADS=off NTL_EXCEPTIONS=off NTL_GMP_LIP=on CXXFLAGS="-I$nemovdir/include" LDFLAGS="-L$nemovdir/lib -Wl,-rpath,$nemovdir/lib"`)
run(`make -j4`)
run(`make install`)
rm(tmp; recursive=true)

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

run(`$srcs/autogen.sh`)

# out of source-tree building:
try 
   mkdir(joinpath(wdir, "Singular_build"))
catch
end

cd(joinpath(wdir, "Singular_build"))
if debug_build
   run(`$srcs/configure --prefix=$vdir --disable-static --disable-p-procs-static --enable-p-procs-dynamic --enable-shared --with-gmp=$nemovdir --with-flint=$nemovdir --with-ntl=$vdir --without-python --with-readline=no --disable-gfanlib`)
else
   run(`$srcs/configure --prefix=$vdir --disable-static --disable-p-procs-static --enable-p-procs-dynamic --enable-shared --with-gmp=$nemovdir --with-flint=$nemovdir --with-ntl=$vdir --without-python --with-readline=no --disable-gfanlib --with-debug --enable-debug --disable-optimizationflags`)
end
run(`make -j4`)
run(`make install`)

cd(wdir)

push!(Libdl.DL_LOAD_PATH, "$pkgdir/local/lib")

cd(oldwdir)
