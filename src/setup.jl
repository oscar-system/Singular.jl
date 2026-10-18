module Setup

import Singular_jll
import lib4ti2_jll
import libsingular_julia_jll

import BinaryWrappers
import CxxWrap
import FileWatching: Pidfile
import Libdl
import Pkg
import Pkg.Artifacts

const lib4ti2_binpath = BinaryWrappers.@generate_wrappers(lib4ti2_jll)

# make sure Singular can find the wrappers
function __init__()
    ENV["PATH"] = lib4ti2_binpath * ":" * ENV["PATH"]
end

#
# regenerate src/libraryfuncdictionary.jl by parsing the list
# of library functions exported by Singular
#
function regenerate_libraryfuncdictionary(prefixpath)

    library_dir = get(ENV, "SINGULAR_LIBRARY_DIR", abspath(prefixpath, "share", "singular", "LIB"))
    filenames = filter(x -> endswith(x, ".lib"), readdir(library_dir))

    #=
      Loops over all libraries and executes libparse on each.
      The first three lines of the libparse output are general information
      about the library, so we ignore it. We are only interested in the
      first column (library name) and the third column (globally exposed or not).
      All other columns (containing info such as line numbers, library name, etc)
      are ignored.
    =#
    dict = Dict{Symbol, Vector{Tuple{String, String}}}()
    cd(abspath(prefixpath, "bin")) do
       for libfile in filenames
           libname = Symbol(libfile[1:end - 4]) # strip the '.lib' suffix
           full_path = joinpath(library_dir, libfile)
           output = read(`$(Singular_jll.libparse()) $full_path`, String)
           libs_splitted = split(output,"\n")[4:end - 1]
           libs_splitted = [split(i, " ", keepempty = false) for i in libs_splitted]
           dict[libname] = [(j[1], j[3]) for j in libs_splitted]
       end
    end
    return dict
end

const libraryfunctiondictionary = regenerate_libraryfuncdictionary(Singular_jll.find_artifact_dir())

function jll_artifact_dir(the_jll::Module)
    artifacts_toml = joinpath(dirname(dirname(Base.pathof(the_jll))), "StdlibArtifacts.toml")

    # If this file exists, it's a stdlib JLL and we must download the artifact ourselves
    if isfile(artifacts_toml)
        # the artifact name is always equal to the module name minus the "_jll" suffix
        name = replace(string(nameof(the_jll)), "_jll" => "")
        meta = Pkg.Artifacts.artifact_meta(name, artifacts_toml)
        hash = Base.SHA1(meta["git-tree-sha1"])
        if !Pkg.Artifacts.artifact_exists(hash)
            dl_info = first(meta["download"])
            Pkg.Artifacts.download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return Pkg.Artifacts.artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by GMP_jll
    return the_jll.find_artifact_dir()
end

# What a cached build in `deps/install-*` must match to be reusable: the bundled
# sources, and the Singular it was compiled against. An overridden Singular_jll
# keeps the same path when it is rebuilt, hence the mtime.
function build_stamp(src_hash::AbstractString, singular_prefix::AbstractString)
   return string(src_hash, "\n", singular_prefix, "\n", mtime(Singular_jll.Singular_path), "\n")
end

function build_code(src_hash)
   depsdir = abspath(joinpath(@__DIR__, "..", "deps"))
   srcdir = joinpath(depsdir, "src")

   # encode the Julia version into into the build & install paths
   # as different Julia versions tend to not be ABI compatible
   builddir = joinpath(depsdir, "build-$VERSION")
   installdir = joinpath(depsdir, "install-$VERSION")

   gmp_prefix = jll_artifact_dir(Singular_jll.GMP_jll)
   singular_prefix = jll_artifact_dir(Singular_jll)
   stamp = build_stamp(src_hash, singular_prefix)

   # check if we already built the code and if so, just use that
   lib_path = joinpath(installdir, "lib", "libsingular_julia.$(Libdl.dlext)")
   stamp_path = joinpath(installdir, "lib", "libsingular_julia.treehash")
   try
      if read(stamp_path, String) == stamp
         @info "Using already compiled bundled C++ code"
         return lib_path
      end
   catch
   end

   # TODO: check if this is a dev version of Singular.jl; if so, continue. If
   # not, i.e. if this is a regular version, then refuse to run


   @info "Compiling libsingular_julia ..."

   try
      run(`cmake --version`)
   catch
      error("Could not locate 'cmake' binary, needed to build libsingular-julia")
   end

   JlCxx_DIR = joinpath(CxxWrap.prefix_path(), "lib", "cmake", "JlCxx")
   julia_exec = joinpath(Sys.BINDIR, Base.julia_exename())

   Pidfile.mkpidlock("$installdir.lock"; stale_age=60) do
      # delete any previous build or install artifacts
      rm(builddir; force=true, recursive=true)
      rm(installdir; force=true, recursive=true)

      withenv("JULIA_LOAD_PATH" => nothing) do
         run(`cmake
             -DJulia_EXECUTABLE=$julia_exec
             -DJlCxx_DIR=$(JlCxx_DIR)
             -Dgmp_prefix=$gmp_prefix
             -DSingular_PREFIX=$(singular_prefix)
             -DCMAKE_INSTALL_PREFIX=$(installdir)
             -DCMAKE_BUILD_TYPE=Release
             -S $srcdir
             -B $builddir
            `)

         run(`cmake
             --build $builddir
             --config Release
             --target install
             --
             -j$(Sys.CPU_THREADS)
            `)
      end

      write(stamp_path, stamp)
   end

   return lib_path

end

"""
    singular_jll_is_overridden()

Return whether Singular_jll resolves to something other than its own artifact,
either via an `Overrides.toml` or via a dev'ed Singular_jll carrying an
`override` directory.
"""
function singular_jll_is_overridden()
   pkgid = Base.PkgId(Singular_jll)

   Pkg.Artifacts.query_override(pkgid.uuid, "Singular") === nothing || return true

   # JLLWrappers resolves a dev'ed JLL to `<package dir>/override` if present
   src = Base.locate_package(pkgid)
   src === nothing && return false
   return isdir(joinpath(dirname(dirname(src)), "override"))
end

function locate_libsingular()
   jll_path = joinpath(libsingular_julia_jll.find_artifact_dir(), "lib", "libsingular_julia.treehash")
   jll_hash = chomp(read(jll_path, String))

   # compare the sources used to build libsingular_julia_jll with bundled copies
   # by comparing tree hashes
   src_hash = bytes2hex(Pkg.GitTools.tree_hash(joinpath(@__DIR__, "..", "deps", "src")))

   # the uuid is for Singular.jl
   pkginfo = get(Pkg.dependencies(), Base.PkgId(parentmodule(Setup)).uuid, nothing)

   if singular_jll_is_overridden()
      # The libsingular_julia in the JLL was compiled against the Singular in
      # Singular_jll. An overridden Singular_jll may differ in ABI -- notably in
      # the omalloc page size -- and a mismatch is not diagnosed, it crashes at
      # runtime. So build against the Singular actually in use.
      @info "Singular_jll is overridden"
      path = build_code(src_hash)
   elseif jll_hash == src_hash || (pkginfo !== nothing && pkginfo.is_tracking_registry)
       # if the tree hashes match then we use the JLL
       # also if we are using a released Singular.jl version
       path = libsingular_julia_jll.get_libsingular_julia_path()
   else
      @info "Bundled C++ sources don't match libsingular_julia_jll"
      path = build_code(src_hash)
   end
   @debug "Use libsingular_julia from $path"
   return path
end

end # module
