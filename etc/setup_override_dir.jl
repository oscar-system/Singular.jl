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
overwrite_allow = false
verbose = false
# debugmode = false
left_ARGS = String[]
while !isempty(ARGS)
   arg = popfirst!(ARGS)
   if arg == "--no-configure"
      global run_configure = false
   elseif arg == "--yes"
      global overwrite_allow = true
   elseif arg == "--verbose"
      global verbose = true
   # elseif arg == "--debug"
   #    global debugmode = true
   else
      push!(left_ARGS, arg)
   end
end


# validate arguments
isdir(singular_prefix) || error("The given Singular prefix '$(singular_prefix)' is not a valid directory")
if ispath(prefix)
   if !overwrite_allow
      print("The given installation prefix '$(prefix)' already exists. Overwrite? [Y/n] ")
      overwrite_allow = lowercase(readline()) in ["y", "yes", ""]
   end
   if overwrite_allow
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
Pkg.add(["JLLPrefixes"])
Pkg.instantiate()

using JLLPrefixes

cd(build_dir)

if run_configure
   @info "Configuring Singular in $(build_dir) for $(prefix)"

   deps = ["GMP_jll", "FLINT_jll", "cddlib_jll", "MPFR_jll"]
   artifact_paths = collect_artifact_paths(deps)
   deps_path = mktempdir(; cleanup=false)
   deploy_artifact_paths(deps_path, artifact_paths)

   extraargs = [
        "CPPFLAGS=-I$(joinpath(deps_path, "include"))",
        "LDFLAGS=-L$(joinpath(deps_path, "lib"))",
   ]

   configure_cmd = `$(singular_prefix)/configure
       --prefix=$(prefix)
       --with-libparse
       --enable-shared
       --disable-static
       --enable-p-procs-static
       --disable-p-procs-dynamic
       --enable-gfanlib
       --without-readline
       --without-ntl
       --with-gmp=$(deps_path)
       --with-flint=$(deps_path)
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
       $(left_ARGS)
       `

   verbose && @show configure_cmd

   # TODO: redirect the output of configure into a log file
   @show run(configure_cmd)
end


@info "Building Singular in $(build_dir)"
run(`make -j$(Sys.CPU_THREADS) $(verbose ? "V=1" : "")`)

@info "Installing Singular to $(prefix)"
run(`make install`)
