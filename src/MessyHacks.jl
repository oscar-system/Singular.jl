#=
   messy hack #0: (not really that messy)

   provide wrappers for working with Singular's global variables
=#

export with_degBound, with_multBound

@doc raw"""
    with_degBound(f, degb::Integer)

Evaluate and return `f()` with the Singular global setting `degBound = degb`. The
value of `degBound` is automatically restored upon return; the effect is only a
local one on `f()`. The value `degBound = 0` corresponds to no degree bound in
Singular and this is the starting value.
"""
function with_degBound(f, degb::Integer)
   old_degb = libSingular.set_degBound(Cint(degb))
   local g = nothing
   try
      g = f()
   finally
      libSingular.set_degBound(old_degb)
   end
   return g
end

@doc raw"""
    with_multBound(f, mu::Integer)

Evaluate and return `f()` with the Singular global setting `multBound = mu`. The
value of `multBound` is automatically restored upon return; the effect is only a
local one on `f()`. The value `multBound = 0` corresponds to no multiplicity
bound in Singular and this is the starting value.
"""
function with_multBound(f, mu::Integer)
   old_mu = libSingular.set_multBound(Cint(mu))
   local g = nothing
   try
      g = f()
   finally
      libSingular.set_multBound(old_mu)
   end
   return g
end

for (name, str) in [(:with_fastHC, "OPT_FASTHC")
                    (:with_infRedTail, "OPT_INFREDTAIL")
                    (:with_lazy, "OPT_OLDSTD")
                    (:with_length, "V_LENGTH")
                    (:with_notBuckets, "OPT_NOT_BUCKETS")
                    (:with_prot, "OPT_PROT")
                    (:with_qringNF, "V_QRING")
                    (:with_redTail, "OPT_REDTAIL")
                    (:with_redThrough, "OPT_REDTHROUGH")]
   @eval begin
      function ($name)(f, flag::Bool)
         old_flag = libSingular.set_option($str, flag)
         local g = nothing
         try
            g = f()
         finally
            libSingular.set_option($str, old_flag)
         end
         return g
      end

      export $name
   end
end

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
      Nemo.ZZRing => Nemo.ZZRingElem,
      Nemo.QQField => Nemo.QQFieldElem,
      Nemo.ZZPolyRing => Nemo.ZZPolyRingElem,
      Nemo.fqPolyRepField => Nemo.fqPolyRepFieldElem,
      Nemo.FqPolyRepField => Nemo.FqPolyRepFieldElem,
      Nemo.AbsSimpleNumField => Nemo.AbsSimpleNumFieldElem
   ],

   one_parameter_types = [
      AbstractAlgebra.Generic.LaurentSeriesRing =>
         AbstractAlgebra.Generic.LaurentSeriesRingElem,
      AbstractAlgebra.Generic.FracField =>
         AbstractAlgebra.Generic.FracFieldElem
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
function all_singular_symbols(R::Nemo.Ring)
   # TODO not correct because AA doesn't have this notion (yet).
   # Doesn't matter to Singular.jl because wrapped coeff rings don't work in
   # the interpreter anyways
   return Symbol[]
end

function all_singular_symbols(R::N_GField)
   return singular_symbols(R)
end

function all_singular_symbols(R::N_AlgExtField)
   return all_singular_symbols(parent(R.minpoly))
end

function all_singular_symbols(R::N_FField)
   return vcat(all_singular_symbols(base_ring(R)), singular_symbols(R))
end

function all_singular_symbols(R::PolyRingUnion)
   return vcat(all_singular_symbols(base_ring(R)), singular_symbols(R))
end

function is_bad_name(x::String)
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
   if is_bad_name(x)
      x = replace(replace(replace(x, "[" => "_"), "]" => ""), "," => "_")
      if is_bad_name(x)
         x = replace(replace(replace(x, "(" => "_"), ")" => ""), "-" => "m")
         if is_bad_name(x)
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
