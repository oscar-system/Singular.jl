# libsingular-julia

This is the C++ library accompanying [Singular.jl](https://github.com/oscar-system/Singular.jl).
It implements the C++ interface from julia to Singular using [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl) and [libcxxwrap-julia](https://github.com/JuliaInterop/libcxxwrap-julia).

This library is compiled into a Julia artifact and shipped to the user via
the Julia JLL package
[libsingular_julia_jll](https://github.com/JuliaBinaryWrappers/libsingular_julia_jll.jl).
This JLL package in turn is generated from the sources in this repository and
the build recipe at
<https://github.com/JuliaPackaging/Yggdrasil/tree/master/L/libsingular_julia>.


## Quick way to test changes in here

To test changes you make to the C++ code here, the quickest way is to use the
`run.jl` script bundled in this repository. As a prerequisite, you need a
working C++ compiler and of course Julia.

Start Julia from this repository as follows:

    julia --project=. run.jl

This will install a few required Julia packages then build all C++ code, and
finally start a Julia session with an artifact override in place which ensures
that libsingular_julia_jll picks up the copy of the C++ code that was just
compiled.

To verify the override works, check `libsingular_julia_jll.artifact_dir`; it
should point at a subdirectory of the current directory.

You can then use the Julia package manager to install or dev `Singular.jl`,
and run its test suite or perform other tests.


## Building

Compiling `libsingular-julia` from source requires a C++ enabled compiler.

The easiest way to build it is to execute the Julia script `build.jl` (which
in turn is also used by `run.jl`). For this you need to execute it in a Julia
environment in which `Singular_jll` and `CxxWrap` are installed.

Alternatively, you can also link it against your own `libcxxwrap-julia`
installation and a Singular installation, but this is more work; an
incantation like the following will do it:

```bash
git clone https://github.com/oscar-system/libsingular-julia \
cmake -DJulia_PREFIX=/home/user/path/to/julia \
      -DSingular_PREFIX=/home/user/path/to/singular
      -DCMAKE_INSTALL_PREFIX=home/user/prefix/for/libsingular-julia \
      -DJlCxx_DIR=/home/user/path/to/libcxxwrap-julia/lib/cmake/JlCxx \
      -DCMAKE_BUILD_TYPE=Release \
      -S libsingular-julia -B build \

cmake --build build --config Release --target install -- -j${nproc}
```


### Overriding the default artifacts for Singular.jl

The `run.jl` script takes care of everything described below, but we document
it in case you need to do any of this manually for some reason.

Put the following into `~/.julia/artifacts/Overrides.toml` to replace the `libsingular-julia` artifact:

```toml
[ae4fbd8f-ecdb-54f8-bbce-35570499b30e]
libsingular_julia = "/home/user/prefix/for/libsingular-julia"
```

If you were using custom versions of `Singular` and/or `libcxxwrap-julia`,
then their directories used during the build need to be added as well, e.g.:

```toml
[bcd08a7b-43d2-5ff7-b6d4-c458787f915c]
Singular = "/home/user/path/to/Singular"

[3eaa8342-bff7-56a5-9981-c04077f7cee7]
libcxxwrap_julia = "/home/user/path/to/libcxxwrap-julia"
```
