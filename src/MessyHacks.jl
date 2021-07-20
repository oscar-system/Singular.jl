#=
   messy hack #1:

   Getting the nemo ring element out of the darn n_unknown struct:

   all types would be handled by the general signature

   function (R::parent_type(T))(a::n_Ring{T}) where {T <: Nemo.RingElem}

   but this is an error, so the following nonsense ensues.

   If you want your stuff to work with the R(a::n_Ring{T}) syntax, you will
   have to add it to the {zero|one}_parameter_types list below or make your
   own special definition.
=#

let zero_parameter_types = [
      Nemo.FlintIntegerRing => Nemo.fmpz,
      Nemo.FlintRationalField => Nemo.fmpq,
      Nemo.FmpzPolyRing => Nemo.fmpz_poly,
      Nemo.FqNmodFiniteField => Nemo.fq_nmod,
      Nemo.FqFiniteField => Nemo.fq,
      Nemo.AnticNumberField => Nemo.nf_elem
   ],

   one_parameter_types = [
      AbstractAlgebra.Generic.LaurentSeriesRing =>
         AbstractAlgebra.Generic.LaurentSeriesRingElem,
      AbstractAlgebra.Generic.FracField =>
         AbstractAlgebra.Generic.Frac
   ]

   for (A, B) in zero_parameter_types
      nu = (A <: Nemo.Field) ? :(Singular.n_FieldElem) : :(Singular.n_RingElem)
      @eval begin
         function (R::($A))(a::($nu){($B)})
            GC.@preserve a begin
               return R(libSingular.julia(libSingular.cast_number_to_void(a.ptr)))
            end
         end
      end
   end

   for (A, B) in one_parameter_types
      nu = (A <: Nemo.Field) ? :(Singular.n_FieldElem) : :(Singular.n_RingElem)
      @eval begin
         function (R::($A){T})(a::($nu){($B){T}}) where T <: RingElement
            GC.@preserve a begin
               return R(libSingular.julia(libSingular.cast_number_to_void(a.ptr)))
            end
         end
      end
   end
end
