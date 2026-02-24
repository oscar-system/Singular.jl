# This Julia script sets up a directory with Singular compiled from a local path,
# then "installs" this GAP for use by the `run_with_override.jl` script


#
# parse arguments
#
length(ARGS) >= 1 || error("must provide path of Singular source directory as first argument")
length(ARGS) >= 2 || error("must provide path of destination directory as second argument")
singular_prefix = popfirst!(ARGS)
prefix = popfirst!(ARGS)
if length(ARGS) > 0 && ARGS[1] == "--debug"
  debugmode = true
  popfirst!(ARGS)
else
  debugmode = false
end

# TODO: should the user be allowed to provide a tmp_gap_build_dir ? that might
# be handy for incremental updates


# validate arguments
isdir(singular_prefix) || error("The given GAP prefix '$(singular_prefix)' is not a valid directory")
if ispath(prefix)
    error("installation prefix '$(prefix)' already exists, please remove it before running this script")
    # TODO: prompt the user for whether to delete the dir or abort
end

# convert into absolute paths
mkpath(prefix)
prefix = abspath(prefix)
singular_prefix = abspath(singular_prefix)

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
# locate GMP headers and the Julia executable for use by the GAP build system
#
deps_and_dirs = Dict(
   GMP_jll => gmp_artifact_dir(),
   FLINT_jll => FLINT_jll.find_artifact_dir(),
   cddlib_jll => cddlib_jll.find_artifact_dir(),
   MPFR_jll => MPFR_jll.find_artifact_dir(),
)

#
# create a temporary directory for the build
#
tmp_gap_build_dir = mktempdir(; cleanup = false)
cd(tmp_gap_build_dir)

#
# configure and build Singular
#
@info "Configuring GAP in $(tmp_gap_build_dir) for $(prefix)"

extraargs = String[]
cppflags = String[]
ldflags = String[]

for (jll, dir) in deps_and_dirs
    push!(cppflags, "-I$(dir)/include")
    push!(ldflags, "-L$(dir)/lib")
end

if !isempty(cppflags)
   push!(extraargs, "CPPFLAGS=$(join(cppflags, " "))")
end
if !isempty(ldflags)
   push!(extraargs, "LDFLAGS=$(join(ldflags, " "))")
end

if debugmode
    # compile GAP in debug mode (enables many additional assertions in the kernel)
    # and disable optimizations, so that debugging the resulting binary with gdb or lldb
    # gets easier
    @info "Debug mode is enabled"
    append!(extraargs, ["CFLAGS=-g", "CXXFLAGS=-g", "--enable-debug"])
end


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


@info "Building Singular in $(tmp_gap_build_dir)"

# complete the build
run(`make -j$(Sys.CPU_THREADS)`)

@info "Installing Singular to $(prefix)"
run(`make install`)
