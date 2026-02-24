#
# parse arguments
#
length(ARGS) >= 1 || error("must provide path of GAP override directory as first argument")
gapoverride = popfirst!(ARGS)

isdir(gapoverride) || error("The given override path '$(gapoverride)' is not a valid directory")
gapoverride = abspath(gapoverride)

#
#
#
@info "Install needed packages"
using Pkg
Pkg.develop(path=dirname(dirname(@__FILE__)))
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

    # we need to make sure that precompilation is run again with the override in place
    # (just running Pkg.precompile() does not seem to suffice)
    run(`touch $(Base.locate_package(pkgid))`)
end

tmpdepot = mktempdir(; cleanup=true)
@info "Created temporary depot at $(tmpdepot)"

# create override file for Singular_jll
add_jll_override(tmpdepot, "Singular", gapoverride)

# prepend our temporary depot to the depot list...
withenv("JULIA_DEPOT_PATH"=>tmpdepot*":"*join(DEPOT_PATH, ":")) do

    # ... and start Julia, by default with the same project environment
    run(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(ARGS)`)
end
