# libsingular-julia

This is the C++ library accompanying [Singular.jl](https://github.com/oscar-system/Singular.jl).
It implements the C++ interface from julia to Singular using [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl) and [libcxxwrap-julia](https://github.com/JuliaInterop/libcxxwrap-julia).

For Singular.jl versions less than 0.3.2 this was included in Singular.jl but version 0.4 will probably use this library as a separate artifact.

## Building

Compiling `libsingular-julia` from source requires a C++ enabled compiler, a `libcxxwrap-julia` installation and a Singular installation.

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

Put the following into `~/.julia/artifacts/Overrides.toml` to replace the `libsingular-julia` artifact:

```toml
[ae4fbd8f-ecdb-54f8-bbce-35570499b30e]
libingular_julia = "/home/user/prefix/for/libsingular-julia"
```

Overrides for `Singular` and `libcxxwrap-julia` with the directories used during the build need to be added as well, e.g.:

```toml
[bcd08a7b-43d2-5ff7-b6d4-c458787f915c]
Singular = "/home/user/path/to/Singular"

[3eaa8342-bff7-56a5-9981-c04077f7cee7]
libcxxwrap_julia = "/home/user/path/to/libcxxwrap-julia"
```
