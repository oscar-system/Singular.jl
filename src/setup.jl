module Setup

using Singular_jll
import lib4ti2_jll
using BinaryWrappers


const lib4ti2_binpath = @generate_wrappers(lib4ti2_jll)

# make sure Singular can find the wrappers
function __init__()
    ENV["PATH"] = lib4ti2_binpath * ":" * ENV["PATH"]
end

#
# regenerate src/libraryfuncdictionary.jl by parsing the list
# of library functions exported by Singular
#
function regenerate_libraryfuncdictionary(prefixpath)

    library_dir = get(ENV, "SINGULAR_LIBRARY_DIR", abspath(prefixpath, "share", "singular", "LIB"))
    filenames = filter(x -> endswith(x, ".lib"), readdir(library_dir))
    output_filename = abspath(@__DIR__, "libraryfuncdictionary.jl")

    @info "Regenerating $(output_filename)"

    #=
      Loops over all libraries and executes libparse on it.
      The first three lines of the libparse output are general information
      about the library, so we ignore it. We are only interested in the
      first column (library name) and the third column (globally exposed or not).
      All other columns (containing info such as line numbers, library name, etc)
      are ignored.
    =#
    io = IOBuffer()
    cd(abspath(prefixpath, "bin")) do
       println(io, "libraryfunctiondictionary = Dict(")
       for libfile in filenames
           full_path = joinpath(library_dir, libfile)
           output = Singular_jll.libparse() do exe
               read(`$exe $full_path`, String)
           end
           libs_splitted = split(output,"\n")[4:end - 1]
           libs_splitted = [split(i, " ", keepempty = false) for i in libs_splitted]
           println(io, ":$(libfile[1:end - 4]) => [")
           for j in libs_splitted
               println(io, """[ "$(j[1])", "$(j[3])" ],""")
           end
           println(io, "],\n")
       end
       println(io, ")\n")
    end
    data = String(take!(io))
    if !isfile(output_filename) || read(output_filename, String) != data
       write(output_filename, data)
    end
end

regenerate_libraryfuncdictionary(Singular_jll.find_artifact_dir())

end # module
