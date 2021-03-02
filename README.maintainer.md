# Directions for updating Singular.jl

(please edit this file as the process changes)

Here is a record of how we solved
oscar-system/Singular.jl#261
by getting the commit
Singular/Singular@48c4ba6
into `Singular.jl`.

`Singular.jl` depends on two binary artifacts:
`Singular_jll` for the actual binaries of Singular, and
`libsingular_julia_jll` for the cxx-wrapped portion.

First, the pull request
JuliaPackaging/Yggdrasil#1622
bumps the SHA to an appropriate version of the Singular sources
and the version `4.1.3p4` in the comment is actually `4.1.3+3` (the 4 is off by one).
After the changes are merged, the repo
https://github.com/JuliaBinaryWrappers/Singular_jll.jl
is automagically updated and `Singular_jll` version `4.1.3+3` is usable from julia.

Next, since the bug fix in Singular is in a header file which also happens
to be included by the cxx wrappers, `libsingular_julia_jll` needs to be updated
as well. This involves making a new release of the repo
https://github.com/oscar-system/libsingular-julia
Since this is now version `0.2.0`, the version and the SHA needs to be updated by
JuliaPackaging/Yggdrasil#1643
After this is merged, the repo
https://github.com/JuliaBinaryWrappers/libsingular_julia_jll.jl
is automagically updated and `libsingular_julia_jll` version `0.2.0` is usable.

Finally, an updated `Singular.jl` is make by bumping the versions in
oscar-system/Singular.jl@56ebe7b
It is important to ping the JuliaRegistrator so that version `0.4.1` is usable.

## new for libsingular_julia process as of 2021:
update libsingular_julia/common.jl and version 1.3
JuliaPackaging/Yggdrasil#2535
once this is merged in the registry (i.e. JuliaRegistries/General#29881), a dummy commit for version 1.4
JuliaPackaging/Yggdrasil#2536
once this is merged in the registry, another dummy commit for version 1.5
JuliaPackaging/Yggdrasil#2537

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
