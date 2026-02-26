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

    # we need to make sure that precompilation is run again with the override in place
    # (just running Pkg.precompile() does not seem to suffice)
    run(`touch $(Base.locate_package(pkgid))`)
end

tmpdepot = mktempdir(; cleanup=true)
@info "Created temporary depot at $(tmpdepot)"

# create override file for Singular_jll
add_jll_override(tmpdepot, "Singular", singularoverride)
run(`touch $(Base.locate_package(Base.identify_package("Singular")))`)

# create a fresh marker file in deps/src so the tree hash changes
libsingular_src_dir = joinpath(dirname(@__DIR__), "deps", "src")
isdir(libsingular_src_dir) || error("Could not find $(libsingular_src_dir)")

marker_prefix = ".recompile-libsingular-julia-"
for filename in readdir(libsingular_src_dir)
    if startswith(filename, marker_prefix)
        rm(joinpath(libsingular_src_dir, filename); force=true)
    end
end

libsingular_src_marker = joinpath(
    libsingular_src_dir,
    marker_prefix * string(time_ns()) * ".txt",
)
write(libsingular_src_marker, "force recompilation marker\n")

# prepend our temporary depot to the depot list...
try
    withenv("JULIA_DEPOT_PATH"=>tmpdepot*":"*join(DEPOT_PATH, ":"), "FORCE_LIBSINGULAR_JULIA_COMPILATION"=>"true") do

        # ... and start Julia, by default with the same project environment
        run(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(ARGS)`)
    end
finally
    rm(libsingular_src_marker; force=true)
end
