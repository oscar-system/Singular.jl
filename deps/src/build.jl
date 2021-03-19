import Singular_jll, CxxWrap, CMake
using Singular_jll.GMP_jll, Pkg.Artifacts

# TODO: use ARGS to specify custom build dir?
builddir = "build"
installdir = abspath("install")
JlCxx_DIR = joinpath(CxxWrap.prefix_path(), "lib", "cmake", "JlCxx")

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
            download_artifact(hash, dl_info["url"], dl_info["sha256"])
        end
        return artifact_path(hash)
    end

    # Otherwise, we can just use the artifact directory given to us by GMP_jll
    return GMP_jll.find_artifact_dir()
end

const gmp_prefix = gmp_artifact_dir()
const singular_prefix = Singular_jll.artifact_dir

rm(builddir; force=true, recursive=true)

run(`$(CMake.cmake)
    -DJulia_EXECUTABLE=$(joinpath(Sys.BINDIR, Base.julia_exename()))
    -Dextra_cppflags=-I$(gmp_prefix)/include
    -Dextra_ldflags=-L$(gmp_prefix)/lib
    -DSingular_PREFIX=$(singular_prefix)
    -DCMAKE_INSTALL_PREFIX=$(installdir)
    -DJlCxx_DIR=$(JlCxx_DIR)
    -DCMAKE_BUILD_TYPE=Release
    -S .
    -B $(builddir)
`)

run(`$(CMake.cmake) --build $(builddir) --config Release --target install -- -j$(Sys.CPU_THREADS)`)
