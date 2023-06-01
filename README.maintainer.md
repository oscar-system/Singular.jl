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
   > ex: <https://github.com/oscar-system/Singular.jl/commit/73af1a6b0f99c11f00837c535db818ca2de7d9a2>

2. After the changes are merged (and before the next `Singular.jl` release), update
   the `libsingular_julia` build script with a new version number and using the
   latest commit SHA for the `master` branch of `Singular.jl`.
   > ex: <https://github.com/JuliaPackaging/Yggdrasil/commit/dd9d775f530b1164dd1cf6135677846dcb25fafc>

3. Wait for this to be merged into Yggdrasil, and then wait for the registry
   to pick up the new version of `libsingular_julia_jll`.
   > ex: <https://github.com/JuliaRegistries/General/commit/00624106c204ccb28fb2c834877fb6005c653ce7>

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

1. Update the Singular build script with the commit SHA of the singular sources (https://github.com/Singular/Singular)
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
