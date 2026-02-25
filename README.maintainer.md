# Directions for updating Singular.jl

`Singular.jl` depends on the Singular kernel and some C++ glue code found in
the `deps/src/` directory. Compiled versions of each are distributed to users
as binary artifacts via the Julia "JLL" packages `Singular_jll` and
`libsingular_julia_jll` respectively.

The build scripts for these JLL packages can be found here:

- <https://github.com/JuliaPackaging/Yggdrasil/blob/master/S/Singular/build_tarballs.jl>
- <https://github.com/JuliaPackaging/Yggdrasil/blob/master/L/libsingular_julia/build_tarballs.jl>

The resulting JLL packages

- <https://github.com/JuliaBinaryWrappers/Singular_jll.jl>
- <https://github.com/JuliaBinaryWrappers/libsingular_julia_jll.jl>

The sources:

- `Singular` sources: <https://github.com/Singular/Singular>
- `libsingular_julia` sources: `deps/src/` directory of `Singular.jl`

## Updating just the C++ wrappers

Suppose just the C++ wrappers need to be updated, without any changes to the
Singular kernel itself.

1. Commit changes to the `deps/src/` directory.
   > ex: <https://github.com/oscar-system/Singular.jl/pull/661>

2. After the changes are merged (and before the next `Singular.jl` release), update
   the `libsingular_julia` build script with a new version number and using the
   latest commit SHA for the `master` branch of `Singular.jl`.
   > ex: <https://github.com/JuliaPackaging/Yggdrasil/pull/6880>

3. Wait for this to be merged into Yggdrasil, and then wait for the registry
   to pick up the new version of `libsingular_julia_jll`.
   > ex: <https://github.com/JuliaRegistries/General/pull/85185>

4. Bump the dependence in `Singular.jl` to whatever version number was used in Step 2.
   > ex: <https://github.com/oscar-system/Singular.jl/commit/f669cea1aa4fc73082d8e8c08a1e33a7c22882ed>

   Version compatibility notation: <https://pkgdocs.julialang.org/v1/compatibility/>

5. Release a new `Singular.jl`. This is done by pinging JuliaRegistrator in the comments of a commit.
   > ex: <https://github.com/oscar-system/Singular.jl/commit/cbc04de26507b386e27d2bd7a23a2f913903d9cc>

After the new version of `Singular.jl` is picked up by the registry, it may be used
in further downstream packages.

## Updating the Singular kernel

Suppose the Singular kernel needs an update. This involves updating both build
scripts because `libsingular_julia_jll` will need to point to the new `Singular_jll`.

1. Update the Singular build script with the commit SHA of the singular sources at <https://github.com/Singular/Singular>
   > ex: <https://github.com/JuliaPackaging/Yggdrasil/commit/1996668f09e89bc8da0cc498c59d905821b150c8>

   In this specific commit, `FLINT_jll` was also updated, but this is not necessary to update singular.
   Any build issues need to be communicated to <https://github.com/Singular/Singular>
   until you get a commit that builds on all targets.

2. Wait for the Yggdrasil merge, and wait for the registry.
   > ex: <https://github.com/JuliaRegistries/General/commit/c73ad5d2ec0483113ffc2eca9b623049a7d8e81a>

3. Update the `libsingular_julia` build scripts with a new version and `Singular_jll` dependency.
   > ex: <https://github.com/JuliaPackaging/Yggdrasil/commit/cc75959bcf85b016625d9e528db277fd4d3b32d9>

At this point, we have a new `libsingular_julia_jll` in the works, and the steps
are essentially Steps 3-5 in the previous section.

4. The usual waiting.
   > ex: <https://github.com/JuliaRegistries/General/commit/5c91a1184ed468d2a014c149f19caf84748fb83f>

5. Bump the `libsingular_julia_jll` and `Singular_jll` dependencies in `Singular.jl`.
   > ex: <https://github.com/oscar-system/Singular.jl/commit/138eeb12b5105e252c086bc45f259704e74c9885>

   This commit also bumps `AbstractAlgebra` and `Nemo`, but this is not necessary
   just to update Singular.

6. Release new `Singular.jl` version.
   > ex: <https://github.com/oscar-system/Singular.jl/commit/12f3b6b8689177d55f3952a8e5f9dccf233be885>

## Updating both `libsingular_julia` and the Singular kernel

Since updating the Singular kernel requires an update to `libsingular_julia`, the
steps here are the same as in the previous section. Just make sure that in
Step 3, the commit SHA used to update the `libsingular_julia` build scripts
contains all of the desired changes to `libsingular_julia`.


## Building a custom `Singular_jll` locally

For testing purposes one may wish to try out `Singular_jll` changes locally before
submitting them as a PR to Yggdrasil. This can be done as shown in the following
shell script:

```shell
# Change into a clone of the Yggdrasil repository
git clone https://github.com/JuliaPackaging/Yggdrasil
cd Yggdrasil

# record the base path
BASEPATH=$(pwd)

# ensure building macOS binaries will work (you can omit this if you only
# want to build for Linux)
export BINARYBUILDER_AUTOMATIC_APPLE=true

# change into the directory containing the `build_tarballs.jl` we want to build
cd S/Singular

# Now `build_tarballs.jl` can be modified, e.g. to pull a different set of
# sources, use different versions of dependencies, etc.

# ensure BinaryBuilder etc. is installed in the right version
# (ideally use the same Julia version as specified in `.ci/Manifest.toml`)
julia --project=$BASEPATH/.ci -e 'using Pkg; Pkg.instantiate()'

# get list of platforms etc.
julia --project=$BASEPATH/.ci build_tarballs.jl --help

# build and deploy the JLL locally. If you omit the comma-separated
# list of PLATFORMS then it will build for *all* platforms
julia --project=$BASEPATH/.ci build_tarballs.jl PLATFORMS --deploy=local
```
