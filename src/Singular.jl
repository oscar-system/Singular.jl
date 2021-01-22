@doc Markdown.doc"""
  Julia package for using the Singular library for commutative and non-commutative algebra, algebraic geometry, and singularity theory.

  For documentation see https://oscar-system.github.io/Singular.jl/latest/index.html

  For more information about Singular see https://www.singular.uni-kl.de/
"""
module Singular

import AbstractAlgebra
using Markdown
using Nemo
using Pkg
using lib4ti2_jll
using Singular_jll
using libsingular_julia_jll

import Base: abs, checkbounds, convert, deepcopy, deepcopy_internal,
             denominator, div, divrem, exponent,
             gcd, gcdx, getindex, hash, inv, isequal, isless, lcm,
             length, mod, numerator, one, reduce, rem, setindex!, show,
             transpose, zero, +, -, *, ==, ^, &, |, <<, >>, ~, <=,
             >=, <, >, //, /, !=

import LinearAlgebra: normalize!, rank

import Statistics: std

using Random: Random, AbstractRNG, SamplerTrivial, SamplerSimple
import Random: rand
using RandomExtensions: RandomExtensions, make, Make2

import AbstractAlgebra: AbstractAlgebra, diagonal_matrix, factor,
                        identity_matrix, kernel, lead, ncols, ngens, nrows, order,
                        preimage, zero_matrix

import Nemo: add!, addeq!, base_ring, canonical_unit,
             change_base_ring, characteristic, check_parent, codomain,
             coeff, coeffs, compose, contains, content, crt,
             deflate, deflation, degree, degrees, derivative, divexact,
             divides, domain, elem_type, evaluate, exponent_vectors, factor_squarefree,
             finish, gcdinv, gen, gens, get_field, intersect, isconstant,
             isgen, ismonomial, inflate, isnegative, isone,
             isterm, isunit, iszero, lift, lc, lt, lm, monomials,
             MPolyBuildCtx, mul!, nvars, ordering, parent_type,
             parent, primpart, promote_rule, push_term!,
             reconstruct, remove, set_field!, sort_terms!,
             symbols, terms, total_degree, valuation,
             var_index, vars, zero!, ResidueRing

export base_ring, elem_type, parent_type, parent

export ResidueRing, PolynomialRing, Ideal, MaximalIdeal, FreeModule

export ZZ, QQ, FiniteField, FunctionField, CoefficientRing, Fp

###############################################################################
#
#   Set up environment / load libraries
#
###############################################################################

const libsingular = Singular_jll.libsingular
const binSingular = Singular_jll.Singular_path

const libsingula_julia = libsingular_julia_jll.libsingular_julia

const libflint = Nemo.libflint
const libantic = Nemo.libantic

mapping_types_reversed = nothing
#
# More helpful error message for users on Windows.
windows_error() = error("""

    This package unfortunately does not run natively under Windows.
    Please install Julia using Windows subsystem for Linux and try again.
    See also https://oscar.computeralgebra.de/install/.
    """)

if Sys.iswindows()
  windows_error()
end

function __init__()
   if Sys.iswindows()
      windows_error()
   end

   # Initialise Singular
   ENV["SINGULAR_EXECUTABLE"] = binSingular
   libSingular.siInit(binSingular)
   # set up Singular parents (we cannot do this before Singular is initialised)

   global ZZ = Integers()
   global QQ = Rationals()

   # done in __init__ since headers must be included first

   global n_Z_2_n_Q = libSingular.n_SetMap(ZZ.ptr, QQ.ptr)
   global n_Q_2_n_Z = libSingular.n_SetMap(QQ.ptr, ZZ.ptr)
   ZZ.refcount += 1
   QQ.refcount += 1

   global ringorder_no = libSingular.ringorder_no
   global ringorder_lp = libSingular.ringorder_lp
   global ringorder_rp = libSingular.ringorder_rp
   global ringorder_dp = libSingular.ringorder_dp
   global ringorder_Dp = libSingular.ringorder_Dp
   global ringorder_ls = libSingular.ringorder_ls
   global ringorder_rs = libSingular.ringorder_rs
   global ringorder_ds = libSingular.ringorder_ds
   global ringorder_Ds = libSingular.ringorder_Ds
   global ringorder_c  = libSingular.ringorder_c
   global ringorder_C  = libSingular.ringorder_C

   global sym2ringorder = Dict{Symbol, libSingular.rRingOrder_t}(
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
   global mapping_types_reversed, casting_functions
   mapping_types_reversed = Dict( i[2] => i[1] for i in libSingular.get_type_mapper() )
   casting_functions = create_casting_functions()

   show_banner = isinteractive() &&
                !any(x->x.name in ["Oscar"], keys(Base.package_locks))

   singular_version_nr=Singular.libSingular.version()
   ver = digits(singular_version_nr, base = 10)
   svn = "$(ver[4]).$(ver[3]).$(ver[2])"
   if ver[1] > 0
     svn *= "p$(ver[1])"
   end
   if show_banner
     println("""Singular.jl, based on
                     SINGULAR                                 /
 A Computer Algebra System for Polynomial Computations       /  Singular.jl: $VERSION_NUMBER
                                                           0<   Singular   : $svn
 by: W. Decker, G.-M. Greuel, G. Pfister, H. Schoenemann     \\
FB Mathematik der Universitaet, D-67653 Kaiserslautern        \\
     """)
   end
end

# pkgdir was added in Julia 1.4
if VERSION < v"1.4"
   pkgdir(m::Core.Module) = abspath(Base.pathof(Base.moduleroot(m)), "..", "..")
end
pkgproject(m::Core.Module) = Pkg.Operations.read_project(Pkg.Types.projectfile_path(pkgdir(m)))
pkgversion(m::Core.Module) = pkgproject(m).version
const VERSION_NUMBER = pkgversion(@__MODULE__)


###############################################################################
#
#   Load Singular Rings/Fields/etc
#
###############################################################################

include("setup.jl")

include("AbstractTypes.jl")

include("LibSingular.jl")

include("Number.jl")

include("Poly.jl")

include("Module.jl")

include("Ideal.jl")

include("Matrix.jl")

include("Vector.jl")

include("Resolution.jl")

include("caller.jl")

include("Meta.jl")

include("Map.jl")

end # module
