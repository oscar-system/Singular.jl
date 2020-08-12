function execute(libparse, path)
  res = libparse() do exe
    cmd = `$exe $path`
    out = Pipe()
    err = Pipe()
  
    process = run(pipeline(ignorestatus(cmd), stdout = out, stderr = err))
    close(out.in)
    close(err.in)
  
    return (stdout = String(read(out)), 
            stderr = String(read(err)),  
            code = process.exitcode)
  end
end

libparsepath = abspath(joinpath(prefixpath, "bin", "libparse"))

if haskey(ENV, "SINGULAR_LIBRARY_DIR")
    library_dir = ENV["SINGULAR_LIBRARY_DIR"]
else
    library_dir = abspath(joinpath(prefixpath, "share", "singular", "LIB"))
end

filenames = filter(x -> endswith(x, ".lib"), readdir(library_dir))

output_filename = abspath(joinpath(@__DIR__, "..", "src", "libraryfuncdictionary.jl"))

#=
  Loops over all libraries and executes libparse on it.
  The first three lines of the libparse output are general information
  about the library, so we ignore it. We are only interested in the
  first column (library name) and the third column (globally exposed or not).
  All other columns (containing info such as line numbers, library name, etc)
  are ignored.
=#
cd(abspath(joinpath(prefixpath, "bin"))) do
    open(output_filename, "w") do outputfile
        println(outputfile, "libraryfunctiondictionary = Dict(")
        for i in filenames
            full_path = joinpath(library_dir, i)
            libs = execute(Singular_jll.libparse, full_path)
            if libs.stderr != ""
                error("from libparse: $(libs.stderr)")
            end
            libs_splitted = split(libs.stdout,"\n")[4:end - 1]
            libs_splitted = [split(i, " ", keepempty = false) for i in libs_splitted]
            println(outputfile, ":$(i[1:end - 4]) => [")
            for j in libs_splitted
                println(outputfile, """[ "$(j[1])", "$(j[3])" ],""")
            end
            println(outputfile, "],\n")
        end
        println(outputfile, ")\n")
    end
end
