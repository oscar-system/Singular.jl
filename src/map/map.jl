export AlgebraHomomorphism, codomain, compose, domain

###############################################################################
#
#   Basic manipulation of Algebra Homomorphisms
#
###############################################################################

function domain(f::salghom)
   return f.domain
end

function codomain(f::salghom)
   return f.codomain
end

###############################################################################
#
#   Basic Operations with Algebra Homomorphisms
#
###############################################################################

function map_ideal(f::salghom, I::sideal)

   if base_ring(I) != f.domain
      error("Ideal is not contained in the domain of the
            algebra homomorphism.")
   end

   return Ideal(f.codomain, libSingular.maMapIdeal(I.ptr, f.domain.ptr, f.image.ptr,
          f.codomain.ptr))
end

function map_poly(f::salghom, p::spoly)

   if parent(p) != f.domain
      error("Polynomial is not contained in the domain of the
            algebra homomorphism.")
   end
   return f.codomain(libSingular.maMapPoly(p.ptr, f.domain.ptr, f.image.ptr,
          f.codomain.ptr))
end

function (f::salghom)(p::spoly)
   return map_poly(f, p)
end

function (f::salghom)(I::sideal)
   return map_ideal(f, I)
end

###############################################################################
#
#   Composition Algebra Homomorphisms
#
###############################################################################

function check_composable(f::salghom, g::salghom)
   codomain(f) != domain(g) && error("Incompatible maps")
end

function compose(f::salghom, g::salghom)
   check_composable(f, g)
   I = g(f.image)
   return AlgebraHomomorphism(f.domain, g.codomain, I)
end

*(f, g) = compose(f, g)


###############################################################################
#
#   Algebra Homomorphism constructor
#
###############################################################################

function AlgebraHomomorphism(D::PolyRing, C::PolyRing, I::sideal) where T <:
         Union{Ring, Field}

   if base_ring(I) != C
      error("Defining equations have to be contained in the codomain.")
   end

   if ngens(I) != nvars(D)
      error("Number of defining equations does not match.")
   end

   if D.base_ring != C.base_ring
      error("Base rings do not match.")
   end

   S = typeof(D.base_ring)
   return salghom{S}(D, C, I)
end

