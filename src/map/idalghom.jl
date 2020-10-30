export codomain, compose, domain, IdentityAlgebraHomomorphism, kernel,
       preimage

###############################################################################
#
#   Basic Operations with Identity Algebra Homomorphisms
#
###############################################################################

function map_ideal(f::Map(SIdAlgHom), I::sideal)

   if base_ring(I) != f.domain
      error("Ideal is not in the domain of the
            algebra homomorphism.")
   end

   return I
end

function map_poly(f::Map(SIdAlgHom), p::spoly)

   if parent(p) != f.domain
      error("Polynomial is not in the domain of the
            algebra homomorphism.")
   end
   return p
end

function (f::SIdAlgHom)(p::spoly)
   return map_poly(f, p)
end

function (f::SIdAlgHom)(I::sideal)
   return map_ideal(f, I)
end

###############################################################################
#
#   I/O for Identity Algebra Homomorphisms
#
###############################################################################

function show(io::IO, M::Map(SIdAlgHom))
   println(io, "Identity Algebra Homomorphism with")
   println(io, "")
   println(io, "Domain: ", domain(M))
   println(io, "")
   println(io, "Defining Equations: ", M.image)
end

###############################################################################
#
#   Composition with Identity Algebra Homomorphisms
#
###############################################################################

function compose(f::Map(SIdAlgHom), g::Map(SAlgHom))
   check_composable(f, g)
   return g
end

function compose(f::Map(SAlgHom), g::Map(SIdAlgHom))
   check_composable(f, g)
   return f
end

function compose(f::Map(SIdAlgHom), g::Map(SIdAlgHom))
   check_composable(f, g)
   return g
end

function check_composable(f::Map(SIdAlgHom), g::Map(SAlgHom))
   codomain(f) != domain(g) && error("Incompatible maps")
end

function check_composable(f::Map(SAlgHom), g::Map(SIdAlgHom))
   codomain(f) != domain(g) && error("Incompatible maps")
end

function check_composable(f::Map(SIdAlgHom), g::Map(SIdAlgHom))
   codomain(f) != domain(g) && error("Incompatible maps")
end

###############################################################################
#
#   Preimage and Kernel for Identity Algebra Homomorphisms
#
###############################################################################

@doc Markdown.doc"""
    preimage(f::AbstractAlgebra.Map(SIdAlgHom), I::sideal)

> Returns the preimage of the ideal $I$ under the identity algebra homomorphism.
"""
function preimage(f::Map(SIdAlgHom), I::sideal)

   if base_ring(I) != f.domain
      error("Ideal is not contained in codomain.")
   end

   return I
end

@doc Markdown.doc"""
   kernel(f::AbstractAlgebra.Map(SIdAlgHom))
> Returns the kernel of the identity algebra homomorphism.
"""
function kernel(f::Map(SIdAlgHom))
   return Ideal(f.domain, )
end

################################################################################
#
#  Identity Algebra Homomorphism Constructor
#
################################################################################

@doc Markdown.doc"""
   Identity_Algebra_Homomorphism(R::PolyRing)
> Constructs the canonical identity algebra homomorphism $id: D --> D$,
> where the $i$-th variable of $D$ is mapped to itself.
"""
IdentityAlgebraHomomorphism(R::PolyRing)  = SIdAlgHom{typeof(R.base_ring)}(R)

domain(f::AbstractAlgebra.Map(Singular.SIdAlgHom)) = f.domain

codomain(f::AbstractAlgebra.Map(Singular.SIdAlgHom)) = f.domain

