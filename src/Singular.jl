module Singular

using Nemo
using Cxx

import Base: abs, deepcopy, den, div, divrem, gcd, gcdx, inv, isequal, isless,
             lcm, mod, num, one, rem, show, zero,
             +, -, *, ==, ^, &, |, $, <<, >>, ~, <=, >=, <, >, //,
             /, !=

import Nemo: add!, addeq!, crt, divexact, elem_type, gcdinv, is_negative,
             isone, iszero, isunit, mul!, needs_parentheses, parent_type,
             parent, reconstruct, show_minus_one, zero!, ResidueRing

export SingularQQ, SingularZZ

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
   cxxinclude(joinpath("coeffs", "rmodulon.h"), isAngled = false)
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

   global const SingularQQ = SingularRationalField()
   global const SingularZZ = SingularIntegerRing()

   global const n_Z_2_n_Q = libSingular.n_SetMap(SingularZZ.ptr, SingularQQ.ptr)
   global const n_Q_2_n_Z = libSingular.n_SetMap(SingularQQ.ptr, SingularZZ.ptr)

   global const ringorder_no = @cxx ringorder_no
   global const ringorder_lp = @cxx ringorder_lp
   global const ringorder_rp = @cxx ringorder_rp
   global const ringorder_dp = @cxx ringorder_dp
   global const ringorder_Dp = @cxx ringorder_Dp
   global const ringorder_ls = @cxx ringorder_ls
   global const ringorder_rs = @cxx ringorder_rs
   global const ringorder_ds = @cxx ringorder_ds
   global const ringorder_Ds = @cxx ringorder_Ds
   global const ringorder_c  = @cxx ringorder_c
   global const ringorder_C  = @cxx ringorder_C

   global const sym2ringorder = Dict{Symbol, Cxx.CppEnum}(
   	  :lex => ringorder_lp,
      :revlex => ringorder_rp, 
   	  :neglex => ringorder_ls,
      :negrevlex => ringorder_rs, 
	  :degrevlex => ringorder_dp,
      :deglex => ringorder_Dp,
	  :negdegrevlex => ringorder_ds,
      :negdeglex => ringorder_Ds,
	  :comp1max => ringorder_c,
      :comp1min => ringorder_C
   )
end

###############################################################################
#
#   Load Singular Rings/Fields/etc
#
###############################################################################

include("AbstractTypes.jl")

include("LibSingular.jl")

include("Number.jl")

include("Poly.jl")

end # module
