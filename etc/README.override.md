# How to use the override scripts

## Using Singular.jl with a different version of Singular than what `Singular_jll` provides

This can be useful for various reasons e.g.,

- you need to test Singular.jl with a newer Singular version, perhaps even its master branch

For this to work, follow these instructions:

1. Obtain a copy of the SINGULAR sources, probably from a clone of the SINGULAR git repository.
   Let's say this is in directory `SINGULARROOT`.

2. Execute `./autoconf.sh` in the SINGULARROOT.

3. Build Singular by executing the `etc/setup_override_dir.jl`
   script. It takes as first argument the SINGULARROOT, and as second argument the places where
   the result shall be installed. I recommend to execute this in a separate
   environment, as it may need to install a few things.

   To give a concrete example you could invoke

        julia --proj=override etc/setup_override_dir.jl $SINGULARROOT /tmp/singular_jll_override

4. Use the `etc/run_with_override.jl` script with the exact same Julia executable
   and the override environment we just prepared.

        julia --proj=override etc/run_with_override.jl /tmp/singular_jll_override

5. This opens a Julia session with the override in effect. You can now e.g. load GAP.jl
   via `using GAP`, or install other packages (such as Oscar) and test with them.
