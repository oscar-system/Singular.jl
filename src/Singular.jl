module Singular

using Nemo
using Cxx

import Base: div, gcd, inv, isequal, isless, lcm, one, show, zero,
             +, -, *, ==, ^, &, |, $, <<, >>, ~, <=, >=, <, >, //,
             /, !=

import Nemo: addeq!, divexact, elem_type, is_negative, isone, iszero, mul!,
             needs_parentheses, parent_type, parent, show_minus_one 

export SingularQQ

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

const pkgdir = realpath(joinpath(dirname(@__FILE__), ".."))
const libsingular = joinpath(pkgdir, "local", "lib", "libSingular")

prefix = joinpath(Pkg.dir("Singular"), "local");

addHeaderDir(joinpath(prefix, "include"), kind = C_System)
addHeaderDir(joinpath(prefix, "include", "singular"), kind = C_System)

function __init__()
   Libdl.dlopen(libsingular, Libdl.RTLD_GLOBAL)

   # include Singular header files

   cxxinclude(joinpath("gmp.h"), isAngled = false)
   cxxinclude(joinpath("omalloc", "omalloc.h"), isAngled = false)
   cxxinclude(joinpath("misc", "intvec.h"), isAngled = false)
   cxxinclude(joinpath("misc", "auxiliary.h"), isAngled = false)
   cxxinclude(joinpath("reporter", "reporter.h"), isAngled = false)
   cxxinclude(joinpath("resources", "feFopen.h"), isAngled = false)
   cxxinclude(joinpath("coeffs", "coeffs.h"), isAngled = false)
   cxxinclude(joinpath("polys", "clapsing.h"), isAngled = false)
   cxxinclude(joinpath("coeffs", "bigintmat.h"), isAngled = false)
   cxxinclude(joinpath("polys", "monomials", "ring.h"), isAngled = false)
   cxxinclude(joinpath("polys", "monomials", "p_polys.h"), isAngled = false)
   cxxinclude(joinpath("polys", "simpleideals.h"), isAngled = false)
   cxxinclude(joinpath("kernel", "GBEngine", "kstd1.h"), isAngled = false) 
   cxxinclude(joinpath("kernel", "GBEngine", "syz.h"), isAngled = false)
   cxxinclude(joinpath("kernel", "ideals.h"), isAngled = false)
   cxxinclude(joinpath("kernel", "polys.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "grammar.h"), isAngled = false) 
   cxxinclude(joinpath("Singular", "libsingular.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "fevoices.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "ipshell.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "ipid.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "subexpr.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "lists.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "idrec.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "tok.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "links", "silink.h"), isAngled = false)
   cxxinclude(joinpath("Singular", "fehelp.h"), isAngled = false)

   # Initialise Singular
   
   binSingular = joinpath(prefix, "bin", "Singular")
   ENV["SINGULAR_EXECUTABLE"] = binSingular

   @cxx siInit(pointer(binSingular))

   # set up Singular parents
   # done in __init__ since headers must be included first

   global SingularQQ = SingularRationalField()
end

###############################################################################
#
#   Load Singular Rings/Fields/etc
#
###############################################################################

include("AbstractTypes.jl")

include("LibSingular.jl")

include("Coeffs.jl")

end # module
