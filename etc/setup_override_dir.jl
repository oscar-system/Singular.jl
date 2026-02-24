# This Julia script sets up a directory with Singular compiled from a local path,
# then "installs" this Singular for use by the `run_with_override.jl` script

#
# parse arguments
#
length(ARGS) >= 1 || error("must provide path of Singular source directory as first argument")
length(ARGS) >= 2 || error("must provide path of destination directory as second argument")
singular_prefix = popfirst!(ARGS)
prefix = popfirst!(ARGS)

if length(ARGS) > 0 && !startswith(ARGS[1], "--")
  build_dir = popfirst!(ARGS)
else
  build_dir = mktempdir(; cleanup = true)
end

run_configure = true
# debugmode = false
for arg in ARGS
   if arg == "--no-configure"
      global run_configure = false
   # elseif arg == "--debug"
   #    global debugmode = true
   else
      error("Unknown argument: $(arg)")
   end
end


# validate arguments
isdir(singular_prefix) || error("The given Singular prefix '$(singular_prefix)' is not a valid directory")
if ispath(prefix)
   print("The given installation prefix '$(prefix)' already exists. Overwrite? [Y/n]")
   answer = readline()
   if lowercase(answer) in ["y", "yes", ""]
      rm(prefix; force=true, recursive=true)
   else
      error("Aborting")
   end
end

# convert into absolute paths
mkpath(prefix)
prefix = abspath(prefix)
singular_prefix = abspath(singular_prefix)
mkpath(build_dir)
build_dir = abspath(build_dir)

#
# Install needed packages
#
@info "Install needed packages"
using Pkg
using Artifacts
Pkg.add(["GMP_jll", "FLINT_jll", "cddlib_jll", "MPFR_jll"])
Pkg.instantiate()

using GMP_jll, FLINT_jll, cddlib_jll, MPFR_jll

# In Julia >= 1.6, there is a "fake" GMP_jll which does not include header files;
# see <https://github.com/JuliaLang/julia/pull/38797#issuecomment-741953480>
function gmp_artifact_dir()
    artifacts_toml = joinpath(dirname(dirname(Base.pathof(GMP_jll))), "StdlibArtifacts.toml")

    # If this file exists, it's a stdlib JLL and we must download the artifact ourselves
    if isfile(artifacts_toml)
        meta = artifact_meta("GMP", artifacts_toml)
        hash = Base.SHA1(meta["git-tree-sha1"])
        if !artifact_exists(hash)
            dl_info = first(meta["download"])
            Pkg.Artifacts.download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by GMP_jll
    return GMP_jll.find_artifact_dir()
end

#
# locate headers for use by the GAP build system
#
deps_and_dirs = Dict(
   GMP_jll => gmp_artifact_dir(),
   FLINT_jll => FLINT_jll.find_artifact_dir(),
   FLINT_jll.OpenBLAS32_jll => FLINT_jll.OpenBLAS32_jll.find_artifact_dir(),
   cddlib_jll => cddlib_jll.find_artifact_dir(),
   MPFR_jll => MPFR_jll.find_artifact_dir(),
)

cd(build_dir)

if run_configure
   #
   # configure and build Singular
   #
   @info "Configuring Singular in $(build_dir) for $(prefix)"

   extraargs = String[]
   cppflags = String[]
   ldflags = String[]

   for (jll, dir) in deps_and_dirs
      push!(cppflags, "-I$(dir)/include")
      push!(ldflags, "-L$(dir)/lib")
   end

   push!(ldflags, "-L"*joinpath(MPFR_jll.find_artifact_dir(), "lib", "julia"))

   push!(ldflags, "-lopenblas")
   push!(ldflags, "-lgfortran")

   if !isempty(cppflags)
      push!(extraargs, "CPPFLAGS=$(join(cppflags, " "))")
   end
   if !isempty(ldflags)
      push!(extraargs, "LDFLAGS=$(join(ldflags, " "))")
   end

   @show(extraargs)

   # TODO: redirect the output of configure into a log file
   @show run(`$(singular_prefix)/configure
       --prefix=$(prefix)
       --with-libparse
       --enable-shared
       --disable-static
       --enable-p-procs-static
       --disable-p-procs-dynamic
       --enable-gfanlib
       --without-readline
       --without-ntl
       --with-gmp=$(deps_and_dirs[GMP_jll])
       --with-flint=$(deps_and_dirs[FLINT_jll])
       --without-python
       --with-builtinmodules=gfanlib,syzextra,customstd,interval,subsets,loctriv,gitfan,freealgebra
       --disable-partialgb-module
       --disable-polymake-module
       --disable-pyobject-module
       --disable-singmathic-module
       --disable-systhreads-module
       --disable-cohomo-module
       --disable-machinelearning-module
       --disable-sispasm-module
       $(extraargs)
       $(ARGS)
       `)
end


@info "Building Singular in $(build_dir)"

# complete the build
run(`make -j$(Sys.CPU_THREADS) V=1`)

@info "Installing Singular to $(prefix)"
run(`make install`)
