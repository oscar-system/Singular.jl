const oldwdir = pwd()
const pkgdir = Pkg.dir("Singular") 

wdir = "$pkgdir/deps"
vdir = "$pkgdir/local"

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
run(`./configure DEF_PREFIX="$vdir" SHARED=on NTL_THREADS=off NTL_EXCEPTIONS=off NTL_GMP_LIP=on CXXFLAGS="-I$vdir/include" LDFLAGS="-L$vdir/lib"`)
run(`make -j4`)
run(`make install`)

# Install Singular

const srcs = joinpath(wdir, "Singular")

# get/update sources
try
  run(`git clone -b spielwiese https://github.com/Singular/Sources.git $srcs`)
catch
  cd(srcs)
  run(`git pull --rebase`)
  cd(joinpath(tmp, ntl, "src"))
end  

run(`$srcs/autogen.sh`)

# out of source-tree building:
cd(mktempdir(tmp))
run(`$srcs/configure --prefix=$vdir --disable-static --disable-p-procs-static --enable-p-procs-dynamic --enable-shared --with-gmp=$vdir --with-flint=$vdir --with-ntl=$vdir --without-python --with-readline=no --disable-gfanlib --with-debug --enable-debug --disable-optimizationflags`)
run(`make -j4`)
run(`make install`)
run(`rm -Rf $tmp`)

cd(wdir)

push!(Libdl.DL_LOAD_PATH, "$pkgdir/local/lib")

cd(oldwdir)
