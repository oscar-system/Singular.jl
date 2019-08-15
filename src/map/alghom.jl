export Algebra_Homomorphism, codomain, compose, domain, kernel, preimage

###############################################################################
#
#   I/O for Algebra Homomorphisms
#
###############################################################################

function show(io::IO, M::Singular.Map(SAlgHom))
   println(io, "Algebra Homomorphism with")
   println(io, "")
   println(io, "Domain: ", domain(M))
   println(io, "")
   println(io, "Codomain: ", codomain(M))
   println(io, "")
   println(io, "Defining Equations: ", gens(M.image))
end

###############################################################################
#
#   Basic Operations with Algebra Homomorphisms
#
###############################################################################

function map_ideal(f::AbstractAlgebra.Map(SAlgHom), I::sideal)

   if base_ring(I) != f.domain
      error("Ideal is not in the domain of the
            algebra homomorphism.")
   end

   return Ideal(f.codomain, libSingular.maMapIdeal(I.ptr, f.domain.ptr,
                f.image.ptr, f.codomain.ptr))
end

function map_poly(f::AbstractAlgebra.Map(SAlgHom), p::spoly)

   if parent(p) != f.domain
      error("Polynomial is not in the domain of the
            algebra homomorphism.")
   end
   return f.codomain(libSingular.maMapPoly(p.ptr, f.domain.ptr, f.image.ptr,
          f.codomain.ptr))
end

function (f::SAlgHom)(p::spoly)
   return map_poly(f, p)
end

function (f::SAlgHom)(I::sideal)
   return map_ideal(f, I)
end

###############################################################################
#
#   Composition Algebra Homomorphisms
#
###############################################################################

@doc Markdown.doc"""
   compose(f::AbstractAlgebra.Map(Singular.SAlgHom),
                         g::AbstractAlgebra.Map(Singular.SAlgHom))
> Returns an algebra homomorphism $h: domain(f) --> codomain(g)$,
> where $h = g(f)$.
"""
function compose(f::AbstractAlgebra.Map(Singular.SAlgHom),
                         g::AbstractAlgebra.Map(Singular.SAlgHom))
   check_composable(f, g)
   I = g(f.image)
   return Algebra_Homomorphism(f.domain, g.codomain, I)
end

###############################################################################
#
#   Preimage and Kernel for Algebra Homomorphisms
#
###############################################################################

@doc Markdown.doc"""
   preimage(f::AbstractAlgebra.Map(SAlgHom), I::sideal)
> Returns the preimage of the ideal $I$ under the algebra homomorphism $f$.
"""
function preimage(f::AbstractAlgebra.Map(SAlgHom), I::sideal)

   if base_ring(I) != f.codomain
      error("Ideal is not contained in codomain.")
   end

   # The following catches an error returned by Singular in iparith.cc
   if (isquotientring(f.domain) && !has_global_ordering(f.domain)) ||
          (isquotientring(f.codomain) && !has_global_ordering(f.codomain))
      error("Algorithm not implemented for local quotient rings.")
   end

   return Ideal(f.domain, Singular.libSingular.maGetPreimage(f.codomain.ptr,
                        f.image.ptr, I.ptr, f.domain.ptr))
end

@doc Markdown.doc"""
   kernel(f::AbstractAlgebra.Map(SAlgHom))
> Returns the kernel of the algebra homomorphism $f$.
"""
function kernel(f::AbstractAlgebra.Map(SAlgHom))
   return preimage(f, Ideal(f.codomain, ))
end
###############################################################################
#
#   Algebra Homomorphism constructor
#
###############################################################################

@doc Markdown.doc"""
   Algebra_Homomorphism(D::PolyRing, C::PolyRing, I::sideal)
> Constructs an algebra homomorphism $f: D --> C$, where the $i$-th variable of
> $D$ is mapped to the $i$-th entry of $I$. $D$ and $C$ must be polynomial
> rings over the same base ring.
"""
function Algebra_Homomorphism(D::PolyRing, C::PolyRing, I::sideal)

   if D.base_ring == Singular.ZZ
      error("Base ring ZZ not implemented.")
   end

   if base_ring(I) != C
      error("Defining equations have to be contained in the codomain.")
   end

   if ngens(I) != nvars(D)
      error("Number of defining equations does not match.")
   end

   if  D.base_ring != C.base_ring
      error("Base rings do not match.")
   end

   return SAlgHom{typeof(D.base_ring)}(D, C, I)
end

domain(f::AbstractAlgebra.Map(Singular.SAlgHom)) = f.domain

codomain(f::AbstractAlgebra.Map(Singular.SAlgHom)) = f.codomain

