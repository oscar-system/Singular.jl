export Algebra_Homomorphism, compose, Identity_Algebra_Homomorphism

###############################################################################
#
#   I/O for Algebra Homomorphisms
#
###############################################################################

function show(io::IO, M::Singular.Map(SAlgHom))
   println(io, "Algebra Homomorphism")
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

function map_ideal(f::Map(SAlgHom), I::sideal)

   if base_ring(I) != f.domain
      error("Ideal is not contained in the domain of the
            algebra homomorphism.")
   end

   return Ideal(f.codomain, libSingular.maMapIdeal(I.ptr, f.domain.ptr,
                f.image.ptr, f.codomain.ptr))
end

function map_poly(f::Singular.Map(SAlgHom), p::spoly)

   if parent(p) != f.domain
      error("Polynomial is not contained in the domain of the
            algebra homomorphism.")
   end
   return f.codomain(libSingular.maMapPoly(p.ptr, f.domain.ptr, f.image.ptr,
          f.codomain.ptr))
end

function (f::Singular.Map(SAlgHom))(p::spoly)
   return map_poly(f, p)
end

function (f::Singular.Map(SAlgHom))(I::sideal)
   return map_ideal(f, I)
end

###############################################################################
#
#   Composition Algebra Homomorphisms
#
###############################################################################

function compose(f::Map(Singular.SAlgHom), g::Map(Singular.SAlgHom))
   check_composable(f, g)
   I = g(f.image)
   return Algebra_Homomorphism(f.domain, g.codomain, I)
end

###############################################################################
#
#   Algebra Homomorphism constructor
#
###############################################################################

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

###############################################################################
#
#   Basic Operations with Identity Algebra Homomorphisms
#
###############################################################################

function map_ideal(f::Map(SIdAlgHom), I::sideal)

   if base_ring(I) != f.domain
      error("Ideal is not contained in the domain of the
            algebra homomorphism.")
   end

   return I
end

function map_poly(f::Singular.Map(SIdAlgHom), p::spoly)

   if parent(p) != f.domain
      error("Polynomial is not contained in the domain of the
            algebra homomorphism.")
   end
   return p
end

function (f::Singular.Map(SIdAlgHom))(p::spoly)
   return map_poly(f, p)
end

function (f::Singular.Map(SIdAlgHom))(I::sideal)
   return map_ideal(f, I)
end
###############################################################################
#
#   I/O for Identity Algebra Homomorphisms
#
###############################################################################

function show(io::IO, M::Map(SIdAlgHom))
   println(io, "Identity Homomorphism with")
   println(io, "")
   println(io, "Domain: ", domain(M))
   println(io, "")
   println(io, "Defining Equations: ", gens(M.image))
end

###############################################################################
#
#   Composition Identity Algebra Homomorphisms
#
###############################################################################

function compose(f::Map(Singular.SIdAlgHom), g::Map(Singular.SAlgHom))
   check_composable(f, g)
   return g
end

function compose(f::Map(Singular.SAlgHom), g::Map(Singular.SIdAlgHom))
   check_composable(f, g)
   return f
end

function compose(f::Map(Singular.SIdAlgHom),
         g::Map(Singular.SIdAlgHom))
   check_composable(f, g)
   return g
end

################################################################################
#
#  Identity Algebra Homomorphism Constructor
#
################################################################################
Identity_Algebra_Homomorphism(R::PolyRing)  = Singular.SIdAlgHom{typeof(R.base_ring)}(R)

codomain(f::Map(Singular.SIdAlgHom)) = domain(f)

