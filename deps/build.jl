using BinaryProvider

using Singular_jll

import CxxWrap
import Libdl
import CMake

using GMP_jll
using MPFR_jll

localprefixpath = joinpath(@__DIR__, "usr")

prefixpath = Singular_jll.artifact_dir

# add shell scripts that startup another julia for some lib4ti2 programs
# singular and libsingular are supposed to at least look in
#  $prefixpath/lib/singular/MOD

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

include("parselibs.jl")
