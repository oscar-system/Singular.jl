export is_quotient_ring, QuotientRing

@doc Markdown.doc"""
    is_quotient_ring(R::PolyRing)

Return `true` if the given ring is the quotient of a polynomial ring with
a non - zero ideal.
"""
function is_quotient_ring(R::PolyRingUnion)
   GC.@preserve R return Bool(Singular.libSingular.rIsQuotientRing(R.ptr))
end

###############################################################################
#
#   Quotient ring construction
#
###############################################################################

function QuotientRing(R::PolyRing{T}, I::sideal{spoly{T}}) where T <: Nemo.RingElem
   R == base_ring(I) || error("parent mismatch")
   I.isGB || error("ideal must be a Groebner basis")
   GC.@preserve R I begin
      ptr = libSingular.make_qring(R.ptr, I.ptr)
      if libSingular.have_error()
         libSingular.rDelete(ptr)
         error("could not construct quotient of $R by $I: "*
                 libSingular.get_and_clear_error())
      end
      S = PolyRing{T}(ptr, base_ring(R), symbols(R))
      return (S, gens(S))
   end
end

function QuotientRing(R::PluralRing{T}, I::sideal{spluralg{T}}) where T <: Nemo.RingElem
   R == base_ring(I) || error("parent mismatch")
   I.isGB || error("ideal must be a Groebner basis")
   I.isTwoSided || error("ideal must be two-sided")
   GC.@preserve R I begin
      ptr = libSingular.make_qring(R.ptr, I.ptr)
      if libSingular.have_error()
         libSingular.rDelete(ptr)
         error("could not construct quotient of $R by $I: "*
                 libSingular.get_and_clear_error())
      end
      S = PluralRing{T}(ptr, base_ring(R), symbols(R))
      return (S, gens(S))
   end
end

function QuotientRing(R::LPRing{T}, I::sideal{slpalg{T}}) where T <: Nemo.RingElem
   R == base_ring(I) || error("parent mismatch")
   I.isGB || error("ideal must be a Groebner basis")
   I.isTwoSided || error("ideal must be two-sided")
   GC.@preserve R I begin
      ptr = libSingular.make_qring(R.ptr, I.ptr)
      if libSingular.have_error()
         libSingular.rDelete(ptr)
         error("could not construct quotient of $R by $I: "*
                 libSingular.get_and_clear_error())
      end
      S = LPRing{T}(ptr, base_ring(R), degree_bound(R), symbols(R))
      return (S, gens(S))
   end
end

