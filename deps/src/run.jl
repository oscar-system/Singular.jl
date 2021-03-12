# This script builds libsingular-julia, then starts a fresh Julia instances
# with an appropriate override set for e.g. running the test suite or whatever
# else you'd like. Arguments are passed on verbatim to the subprocess julia,
# so you can invokes this kinda like you'd invoke Julia

@info "Install needed packages"
using Pkg
Pkg.add(["Singular_jll", "CxxWrap", "CMake", "libsingular_julia_jll"])

@info "Build libsingular-julia"
include("build.jl")


mktempdir() do tmpdepot
    @info "Created temporary depot at $(tmpdepot)"

    # create override file for libsingular_julia_jll
    uuid = string(Base.identify_package("libsingular_julia_jll").uuid)
    mkpath(joinpath(tmpdepot, "artifacts"))
    open(joinpath(tmpdepot, "artifacts", "Overrides.toml"), "w") do f
        write(f, """
        [$(uuid)]
        libsingular_julia = "$(installdir)"
        """)
    end

    # prepend our temporary depot to the depot list...
    withenv("JULIA_DEPOT_PATH"=>tmpdepot*":") do
        # ... and start Julia, by default with the same project environment
        run(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(ARGS)`)

        # TODO: perform some additional steps here, e.g. perhaps
        # verify that `libsingular_julia_jll.artifact_dir` is set right?
    end
end
