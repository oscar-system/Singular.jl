export AlgebraHomomorphism, codomain, compose, domain, kernel, preimage

###############################################################################
#
#   I/O for Algebra Homomorphisms
#
###############################################################################

function AbstractAlgebra.show_map_data(io::IO, M::Map(SAlgHom))
   println(io)
   print(io, "Defining equations: ", M.image)
end

function Base.show(io::IO, f::Map(SAlgHom))
  io = pretty(io)
  if is_terse(io)
    print(io, "Algebra homomorphism")
  else
    print(io, "Hom: ")
    print(terse(io), Lowercase(), domain(f), " -> ")
    print(terse(io), Lowercase(), codomain(f))
  end
end

###############################################################################
#
#   Basic Operations with Algebra Homomorphisms
#
###############################################################################

function map_ideal(f::Map(SAlgHom), I::sideal)

   if base_ring(I) != f.domain
      error("Ideal is not in the domain of the
            algebra homomorphism.")
   end

   ptr = GC.@preserve I f libSingular.maMapIdeal(I.ptr, f.domain.ptr,
                f.ptr, f.codomain.ptr, libSingular.ndCopyMap())
   
   J = Ideal(f.codomain,ptr)
   # compare the orderings of f.domain.ptr and f.codomain.ptr
   domain_ord=Cint[]
   libSingular.rOrdering_helper(domain_ord, f.domain.ptr)
   codomain_ord=Cint[]
   libSingular.rOrdering_helper(codomain_ord, f.codomain.ptr)
   equal_ordering=true
   for i in 1:size(domain_ord,1)
     if domain_ord[i] != codomain_ord[i]
       equal_ordering=false
       break
     end
   end
   if equal_ordering
     J.isGB = I.isGB
   end
   return J
end

function map_poly(f::Map(SAlgHom), p::spoly)

   if parent(p) != f.domain
      error("Polynomial is not in the domain of the
            algebra homomorphism.")
   end
   GC.@preserve p f return f.codomain(libSingular.maMapPoly(p.ptr, f.domain.ptr, f.ptr,
          f.codomain.ptr, libSingular.ndCopyMap()))
end

function (f::SAlgHom)(p::spoly)
   return map_poly(f, p)
end

function (f::SAlgHom)(I::sideal)
   return map_ideal(f, I)
end

###############################################################################
#
#   Composition of Algebra Homomorphisms
#
###############################################################################

@doc raw"""
    compose(f::AbstractAlgebra.Map(Singular.SAlgHom),
                         g::AbstractAlgebra.Map(Singular.SAlgHom))

Return an algebra homomorphism $h: domain(f) \to codomain(g)$,
where $h = g(f)$.
"""
function compose(f::Map(SAlgHom), g::Map(SAlgHom))
   check_composable(f, g)

   n = length(f.image)
   R = g.domain
   S = g.codomain

   GC.@preserve f R g S ptr = libSingular.maMapIdeal(f.ptr, R.ptr, g.ptr, S.ptr,
                                libSingular.ndCopyMap())

   V = Vector{spoly}()

   for i in 1:n
      p = libSingular.getindex(ptr, Cint(i - 1))
      GC.@preserve S push!(V, S(libSingular.p_Copy(p, S.ptr)))
   end

   return AlgebraHomomorphism(f.domain, g.codomain, V)
end

function check_composable(f::Map(SAlgHom), g::Map(SAlgHom))
   codomain(f) != domain(g) && error("Incompatible maps")
end

###############################################################################
#
#   Preimage and Kernel for Algebra Homomorphisms
#
###############################################################################

@doc raw"""
    preimage(f::AbstractAlgebra.Map(SAlgHom), I::sideal)

Return the preimage of the ideal $I$ under the algebra homomorphism $f$.
"""
function preimage(f::Map(SAlgHom), I::sideal)

   if base_ring(I) != f.codomain
      error("Ideal is not contained in codomain.")
   end

   # The following catches an error returned by Singular in iparith.cc
   if (is_quotient_ring(f.domain) && !has_global_ordering(f.domain)) ||
          (is_quotient_ring(f.codomain) && !has_global_ordering(f.codomain))
      error("Algorithm not implemented for local quotient rings.")
   end

   GC.@preserve f I return Ideal(f.domain, Singular.libSingular.maGetPreimage(f.codomain.ptr,
                        f.ptr, I.ptr, f.domain.ptr))
end

@doc raw"""
    kernel(f::AbstractAlgebra.Map(SAlgHom))

Return the kernel of the algebra homomorphism $f$.
"""
function kernel(f::Map(SAlgHom))
   return preimage(f, Ideal(f.codomain, ))
end

###############################################################################
#
#   Algebra Homomorphism constructor
#
###############################################################################

@doc raw"""
    AlgebraHomomorphism(D::PolyRing, C::PolyRing, V::Vector)

Constructs an algebra homomorphism $f: D \to C$, where the $i$-th variable of
$D$ is mapped to the $i$-th entry of $V$. $D$ and $C$ must be polynomial
rings over the same base ring.
"""
function AlgebraHomomorphism(D::PolyRing, C::PolyRing, V::Vector)

   n = length(V)

   if n != nvars(D)
      error("Number of defining equations does not match.")
   end

   if  D.base_ring != C.base_ring
      error("Base rings do not match.")
   end

   return SAlgHom{typeof(D.base_ring)}(D, C, V)
end

domain(f::AbstractAlgebra.Map(Singular.SAlgHom)) = f.domain

codomain(f::AbstractAlgebra.Map(Singular.SAlgHom)) = f.codomain

