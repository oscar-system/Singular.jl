#
# parse arguments
#
length(ARGS) >= 1 || error("must provide path of Singular override directory as first argument")
singularoverride = popfirst!(ARGS)

isdir(singularoverride) || error("The given override path '$(singularoverride)' is not a valid directory")
singularoverride = abspath(singularoverride)

#
#
#
@info "Install needed packages"
using Pkg
Pkg.develop(path=dirname(@__DIR__))
Pkg.add(["Singular_jll"])
Pkg.instantiate()

#
#
#
function add_jll_override(depot, pkgname, newdir)
    pkgid = Base.identify_package("$(pkgname)_jll")
    pkguuid = string(pkgid.uuid)
    mkpath(joinpath(depot, "artifacts"))
    open(joinpath(depot, "artifacts", "Overrides.toml"), "a") do f
        write(f, """
        [$(pkguuid)]
        $(pkgname) = "$(newdir)"
        """)
    end
end

tmpdepot = mktempdir(; cleanup=true)
@info "Created temporary depot at $(tmpdepot)"

# create override file for Singular_jll
add_jll_override(tmpdepot, "Singular", singularoverride)

singular_libdir = joinpath(singularoverride, "lib")
dyld_fallback = let existing = get(ENV, "DYLD_FALLBACK_LIBRARY_PATH", "")
    isempty(existing) ? singular_libdir : existing * ":" * singular_libdir
end

# Use the temporary depot alone: a trailing separator appends only the system
# depots, not ~/.julia. Nothing precompiled against the unoverridden Singular_jll
# is visible, so everything below is built with the override in place. This
# matters because an artifact override by itself does not invalidate a package
# image, and Singular.jl bakes plenty into one: the libsingular_julia path, the
# Singular binary path, the library function dictionary, and the CxxWrap
# wrappers generated from libsingular_julia.
withenv(
    "JULIA_DEPOT_PATH"=>tmpdepot*":",
    "DYLD_FALLBACK_LIBRARY_PATH"=>dyld_fallback,
) do

    # ... make sure all dependencies are installed ...
    run(`$(Base.julia_cmd()) --project=$(Base.active_project()) -e "using Pkg; Pkg.instantiate()"`)
    # ... and start Julia, by default with the same project environment
    run(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(ARGS)`)
end
