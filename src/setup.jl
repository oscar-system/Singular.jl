module Setup

using Singular_jll

prefixpath = Singular_jll.artifact_dir

# add shell scripts that startup another julia for some lib4ti2 programs
# singular and libsingular are supposed to at least look in
#  $prefixpath/lib/singular/MOD

# FIXME: This is a HACK as we modify the JLL directory; better would be to
# create a custom location for these files (e.g. via a MutableArtifact or via
# Scratch.jl), then add that to the PATH env var so that Singular can find
# them

mkpath("$prefixpath/lib/singular/MOD")

write("$prefixpath/lib/singular/MOD/graver", """
#!/bin/sh
julia --startup-file=no -O0 -e 'import lib4ti2_jll
lib4ti2_jll.zsolve() do x
  p=run(ignorestatus(`\$x -G \$ARGS`))
  exit(p.exitcode)
end' -- "\$@"
""")

write("$prefixpath/lib/singular/MOD/hilbert", """
#!/bin/sh
julia --startup-file=no -O0 -e 'import lib4ti2_jll
lib4ti2_jll.zsolve() do x
  p=run(ignorestatus(`\$x -H \$ARGS`))
  exit(p.exitcode)
end' -- "\$@"
""")

write("$prefixpath/lib/singular/MOD/markov", """
#!/bin/sh
julia --startup-file=no -O0 -e 'import lib4ti2_jll
lib4ti2_jll.exe4ti2gmp() do x
  p=run(ignorestatus(`\$x markov \$ARGS`))
  exit(p.exitcode)
end' -- "\$@"
""")

chmod("$prefixpath/lib/singular/MOD/graver", 0o777)
chmod("$prefixpath/lib/singular/MOD/hilbert", 0o777)
chmod("$prefixpath/lib/singular/MOD/markov", 0o777)



"""
    execute(cmd)

Run a command object. Returns a named tuple containing the stdout and stderr
as well as the exit code of the command.
"""
function execute(cmd::Base.Cmd)
    out = IOBuffer()
    err = IOBuffer()
    process = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))
    return (stdout = String(take!(out)),
            stderr = String(take!(err)),
            code = process.exitcode)
end

#
# regenerate src/libraryfuncdictionary.jl by parsing the list
# of library functions exported by Singular
#
function regenerate_libraryfuncdictionary(prefixpath)

    library_dir = get(ENV, "SINGULAR_LIBRARY_DIR", abspath(joinpath(prefixpath, "share", "singular", "LIB")))
    filenames = filter(x -> endswith(x, ".lib"), readdir(library_dir))
    output_filename = abspath(joinpath(@__DIR__, "libraryfuncdictionary.jl"))

    @info "Regenerating $(output_filename)"

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
                libs = Singular_jll.libparse() do exe
                    execute(`$exe $full_path`)
                end
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
end

regenerate_libraryfuncdictionary(Singular_jll.find_artifact_dir())

end # module
