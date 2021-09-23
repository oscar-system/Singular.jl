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


#=
   messy hack #2:

   The names used for singular rings have to work with the singular interpreter,
   which constantly prints and parses objects.
=#

# return an array of all symbols used for printing objects in the ring
function all_symbols(R::Nemo.Ring)
   return Symbol[]
end

function all_symbols(R::N_GField)
   return Symbol[R.S]
end

function all_symbols(R::N_AlgExtField)
   return all_symbols(parent(R.minpoly))
end

function all_symbols(R::Union{N_FField, PolyRing})
   return vcat(all_symbols(base_ring(R)), singular_symbols(R))
end

function isbad_name(x::String)
   if !occursin(r"\A[@a-zA-Z\']([@a-zA-Z\']*[0-9]*_*)*\Z", x)
      return true
   end
   # insert list of reserved identifier names here
   return x == "minpoly"
end

# return an appropriate Symbol version of x that avoids
# bad singular interpreter names and does not conflict with prev.
# gensym generates truly awful names, so cannot use that
function rename_symbol(prev::Base.Set{String}, x::String, def::String)
   if isbad_name(x)
      x = replace(replace(replace(x, "[" => "_"), "]" => ""), "," => "_")
      if isbad_name(x)
         x = replace(replace(replace(x, "(" => "_"), ")" => ""), "-" => "m")
         if isbad_name(x)
            x = def
         end
      end
   end

   x in prev || return Symbol(x)

   i = 0
   while (i += 1) < 10000
      xi = x*"@"*string(i)
      xi in prev || return Symbol(xi)
   end

   i = BigInt(i)
   while true
      xi = x*"@"*string(i)
      xi in prev || return Symbol(xi)
   end
end

function rename_symbol(x::String, def::String)
   return rename_symbol(Base.Set{String}(), x, def)
end

function rename_symbols(prev::Vector{Symbol}, xcands::Vector{String}, def::String)
   prev = Base.Set(map(String, prev))
   realxs = Symbol[]
   for i in xcands
      realx = rename_symbol(prev, i, def)
      push!(realxs, realx)
      push!(prev, String(realx))
   end
   return realxs
end
