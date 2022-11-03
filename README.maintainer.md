# Directions for updating Singular.jl

`Singular.jl` depends on the Singular kernel and various C++ wrappers maintained
in `libsingular_julia`. These come as binary artifacts `Singular_jll` and
`libsingular_julia_jll` respectively.

The build scripts for each of these:

- <https://github.com/JuliaPackaging/Yggdrasil/blob/master/S/Singular/build_tarballs.jl>
- <https://github.com/JuliaPackaging/Yggdrasil/blob/master/L/libsingular_julia/build_tarballs.jl>

The sources:

- <https://github.com/Singular/Singular>
- <https://github.com/oscar-system/libsingular-julia>

## updating libsingular_julia

Suppose just the C++ wrappers need to be updated, without any changes to the
Singular kernel itself.

1. Commit changes to the `libsingular_julia` source
    ex: https://github.com/oscar-system/libsingular-julia/commit/65a12e7fb8546851d3296244aeacc5de1d97af2c

2. Update the `libsingular_julia` build script with a new version number the commit SHA
    ex: https://github.com/JuliaPackaging/Yggdrasil/commit/dd9d775f530b1164dd1cf6135677846dcb25fafc

3. Wait for this to be merged into Yggdrasil, and then wait for the registry to pick up the new version of `libsingular_julia_jll`
    ex: https://github.com/JuliaRegistries/General/commit/00624106c204ccb28fb2c834877fb6005c653ce7

4. Bump the dependence in `Singular.jl` to whatever version number was used in Step 2.
    ex: https://github.com/oscar-system/Singular.jl/commit/f669cea1aa4fc73082d8e8c08a1e33a7c22882ed

   Version compatibility notation: https://pkgdocs.julialang.org/v1/compatibility/

5. Release a new `Singular.jl`. This is done by pinging JuliaRegistrator in the comments of a commit.
    ex: https://github.com/oscar-system/Singular.jl/commit/cbc04de26507b386e27d2bd7a23a2f913903d9cc

After the new version of `Singular.jl` is picked up by the registry, it may be used
in further downstream packages.

Note: if you would like test simultaneous changes to both `libsingular_julia` and
`Singular.jl`, create pull requests against both with the same branch name.
The CI should show some "matching: " entries.
    ex: https://github.com/oscar-system/libsingular-julia/pull/55

## updating the Singular kernel

Suppose the Singular kernel needs an update. This involves updating both build
scripts because `libsingular_julia_jll` will need to point to the new `Singular_jll`.

1. Update the Singular build script with the commit SHA of the singular sources (https://github.com/Singular/Singular)
    ex: https://github.com/JuliaPackaging/Yggdrasil/commit/1996668f09e89bc8da0cc498c59d905821b150c8

   In this specific commit, `FLINT_jll` was also updated, but this is not necessary to update singular.
   Any build issues need to be communicated to https://github.com/Singular/Singular
   until you get a commit that builds on all targets.

2. Wait for the Yggdrasil merge, and wait for the registry.
    ex: https://github.com/JuliaRegistries/General/commit/c73ad5d2ec0483113ffc2eca9b623049a7d8e81a

3. Update the `libsingular_julia` build scripts with a new version and `Singular_jll` dependency.
    ex: https://github.com/JuliaPackaging/Yggdrasil/commit/cc75959bcf85b016625d9e528db277fd4d3b32d9

   Supposedly, Yggdrasil does not like new versions with the same commit SHA.
   Therefore, if `https://github.com/oscar-system/libsingular-julia.git` has not
   changed since the last version, add a dummy commit and use the new SHA.

At this point, we have a new `libsingular_julia_jll` in the works, and the steps
are essentially Steps 3-5 in the previous section.

4. The usual waiting.
    ex: https://github.com/JuliaRegistries/General/commit/5c91a1184ed468d2a014c149f19caf84748fb83f

5. Bump the `libsingular_julia_jll` and `Singular_jll` dependencies in `Singular.jl`.
    ex: https://github.com/oscar-system/Singular.jl/commit/138eeb12b5105e252c086bc45f259704e74c9885

   This commit also bumps `AbstractAlgebra` and `Nemo`, but this is not necessary
   just to update Singular.

6. Release new `Singular.jl` version.
    ex: https://github.com/oscar-system/Singular.jl/commit/12f3b6b8689177d55f3952a8e5f9dccf233be885

## updating both libsingular_julia and the Singular kernel

Since updating the Singular kernel requires an update to `libsingular_julia`, the
steps here are the same as in the previous section. Just make sure that in
Step 3, the commit SHA used to update the `libsingular_julia` build scripts
contains all of the desired changes to `libsingular_julia`.

# Overriding the libsingular_julia_jll artifact

Suppose you have cloned
https://github.com/oscar-system/libsingular-julia
into, say `/my/path`, and would like to test your local changes.

First, you must have a file `.julia/artifacts/Overrides.toml` with the entry

```
[ae4fbd8f-ecdb-54f8-bbce-35570499b30e]
libsingular_julia = "/my/path/libsingular-julia/build"
```

Next, from `/my/path` run `julia build.jl` which will build the package with
your changes.

Finally, the updated `libsingular_julia_jll` should be available the next time
julia starts. Packages that depend on it may need to be precompiled again to
pick up the changes.
