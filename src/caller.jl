# The return of a singular library procedure is always a tuple of normal
# singular values. This tuple can be of length one, indicating one return value.
# Normal singular values (singular lists, polys, ideals, ints, ...) do not have
# tuples in them
#
# Normal singular values are returned from libsingular to julia as a Vector{Any}:
#     [false, ptr, type]
# where ptr is a ptr to the normal singular value and type is int type code.
#
# Tuples of length n > 1 are returned from libsingular as a Vector{Any}:
#     [true, obj1, obj2, ..., objn]
# where each obji is a normal singular value. So, the case of two return values
# from singular (a return of a tuple of length 2) is:
#     [true, [false, ptr1, type1], [false, ptr2, type2]]

# Now, the fun begins when considering how singular lists are returned from
# libsingular. A singular list (which is a normal singular value) of length n
# is returned as Vector{Any}:
#     [obj1, obj2, ..., objn]
# Since true or false is never the representation of a normal singular value,
# these julia objects can be distinguished from the normal (non-list) values
# and tuples of normal values by seeing if the first argument is a Bool.

# Finally the really fun part: When converting everything back to something for
# julia, the following both produce the same julia result [a, b]:
#     - (a,b): a tuple of length 2 returned by a procedure
#     - list(a,b): a list of length 2 returned by a procedure
# Therefore, the user has to know how to interpret the result.

const casting_functions = Dict{Int64, Function}()

function create_casting_functions()
    return Dict(
        mapping_types_reversed[:NUMBER_CMD] =>
            function (vptr, R)
                cast = libSingular.NUMBER_CMD_CASTER(vptr)
                return R.base_ring(cast)
            end
        ,
        mapping_types_reversed[:RING_CMD] =>
            function (vptr, R)
                cast = libSingular.RING_CMD_CASTER(vptr)
                new_ring = create_ring_from_singular_ring(cast)
                return [new_ring, convert_ring_content(libSingular.get_ring_content(cast), new_ring)]
            end
        ,
        mapping_types_reversed[:POLY_CMD] =>
            function (vptr, R)
                cast = libSingular.POLY_CMD_CASTER(vptr)
                return spoly{elem_type(base_ring(R))}(R, cast)
            end
        ,
        mapping_types_reversed[:IDEAL_CMD] =>
            function (vptr, R)
                cast = libSingular.IDEAL_CMD_CASTER(vptr)
                return sideal{elem_type(R)}(R, cast)
            end
        ,
        mapping_types_reversed[:MODUL_CMD] =>
            function (vptr, R)
                cast = libSingular.IDEAL_CMD_CASTER(vptr)
                return smodule{elem_type(R)}(R, cast)
            end
        ,
        mapping_types_reversed[:VECTOR_CMD] =>
            function (vptr, R)
                cast = libSingular.POLY_CMD_CASTER(vptr)
                return svector{elem_type(R)}(R, 1, cast)
            end
        ,
        mapping_types_reversed[:MATRIX_CMD] =>
            function (vptr, R)
                cast = libSingular.MATRIX_CMD_CASTER(vptr)
                return smatrix{elem_type(R)}(R, cast)
            end
        ,
        mapping_types_reversed[:INT_CMD] =>
            function (vptr, R)
                cast = libSingular.INT_CMD_CASTER(vptr)
                return cast::Int
            end
        ,
        mapping_types_reversed[:STRING_CMD] =>
            function (vptr, R)
                cast = libSingular.STRING_CMD_CASTER(vptr)
                return String(cast)
            end
        ,
        mapping_types_reversed[:LIST_CMD] =>
            function (vptr, R)
                cast = libSingular.LIST_CMD_TRAVERSAL(vptr)
                return convert_return(cast, R)
            end
        ,
        mapping_types_reversed[:INTVEC_CMD] =>
            function (vptr, R)
                cast = libSingular.INTVEC_CMD_CASTER(vptr)
                return cast::Vector{Int}
            end
        ,
        mapping_types_reversed[:INTMAT_CMD] =>
            function (vptr, R)
                cast = libSingular.INTMAT_CMD_CASTER(vptr)
                return cast::Matrix{Int}
            end
        ,
        mapping_types_reversed[:BIGINT_CMD] =>
            function (vptr, R)
                cast = libSingular.NUMBER_CMD_CASTER(vptr)
                return libSingular.n_GetMPZ(cast, libSingular.get_coeffs_BIGINT())
            end
        ,
        mapping_types_reversed[:BIGINTMAT_CMD] =>
            function (vptr, R)
                cast = libSingular.BIGINTMAT_CMD_CASTER(vptr)
                return sbigintmat(cast)
            end
        ,
        mapping_types_reversed[:MAP_CMD] =>
            function (vptr, R)
                @warn "returning a map from a Singular library procedure as an ideal" maxlog=1
                # punning in libpolys/polys/simpleideals.h: clear the preimage
                # string of the map and replace it with the rank 1 of the ideal
                libSingular.omFree(unsafe_load(reinterpret(Ptr{Ptr{UInt8}}, vptr), 2))
                unsafe_store!(reinterpret(Ptr{Int}, vptr), 1, 2)
                cast = libSingular.IDEAL_CMD_CASTER(vptr)
                return sideal{elem_type(R)}(R, cast)
            end
        ,
        mapping_types_reversed[:RESOLUTION_CMD] =>
            function (vptr, R)
                cast = libSingular.RESOLUTION_CMD_CASTER(vptr)
                return R(cast, Val(:resolution))  # eh
            end
        )
end

# for the translation of a non-tuple return
# list are possible here, but they are still wrapped in the opaque valueptr
function convert_normal_value(valueptr, typ, R)
    mapper = get(casting_functions, typ) do
       error("unrecognized object with Singular type number $typ\n"*
             "Note that if Singular.jl cannot interpret the type, "*
             "it is doubtful that a interpreter procedure returning "*
             "such a type can be useful to Singular.jl.\n"*
             "Note also that Singular.jl does not support the "*
             "attributes used by the interpreter.\n"*
             "Finally, Singular.lookup_library_symbol can be "*
             "used to fetch the current value of global variables "*
             "stored in the interpreter.")
    end
    return mapper(valueptr, R)
end

function convert_ring_content(value_list, R)
    return Dict(i[2] => convert_normal_value(i[3], i[1], R) for i in value_list)
end

# take ownership of the pointer - not for general users
function create_ring_from_singular_ring(r::libSingular.ring_ptr)
   c = libSingular.rCoeffPtr(r)
   if libSingular.nCoeff_is_Q(c)
      basering = QQ
      T = n_Q
   elseif libSingular.nCoeff_is_Zp(c)
      p = Int(libSingular.n_GetChar(c))
      basering = N_ZpField(p)
      T = n_Zp
   elseif libSingular.nCoeff_is_Z(c)
      basering = ZZ
      T = n_Z
   elseif libSingular.nCoeff_is_GF(c)
      p = Int(libSingular.n_GetChar(c))
      q = Int(libSingular.nfCharQ(c))
      d = round(Int, log(p, q))
      s = Symbol(libSingular.n_ParameterName(0, c))
      basering = N_GField(p, d, s)
      T = n_GF
   elseif libSingular.nCoeff_is_transExt(c)
      p = Int(libSingular.n_GetChar(c))
      npars = libSingular.n_NumberOfParameters(c)
      S = [Symbol(libSingular.n_ParameterName(i-1, c)) for i in 1:npars]
      basering = N_FField(iszero(p) ? QQ : N_ZpField(p), S)
      T = n_transExt
   elseif libSingular.nCoeff_is_algExt(c)
      # first create the univariate transcendental extension
      p = Int(libSingular.n_GetChar(c))
      @assert libSingular.n_NumberOfParameters(c) == 1
      S = [Symbol(libSingular.n_ParameterName(0, c))]
      F = N_FField(iszero(p) ? QQ : N_ZpField(p), S)
      # now create the extension
      minpoly = F(libSingular.algExt_GetMinpoly(c, F.ptr))
      basering = N_AlgExtField(libSingular.nCopyCoeff(c), minpoly)
      T = n_algExt
   elseif libSingular.nCoeff_is_Nemo_Field(c) || libSingular.nCoeff_is_Nemo_Ring(c)
      cf = libSingular.nCopyCoeff(c)
      data_ptr = libSingular.nGetCoeffData(cf)
      R = unsafe_pointer_to_objref(data_ptr)
      basering = CoefficientRing(R)  # FIXME: should we set cache=false ?
      T = elem_type(basering)
   else
      basering = N_UnknownSingularCoefficientRing(libSingular.nCopyCoeff(c))
      T = n_unknownsingularcoefficient
   end
   if libSingular.rIsPluralRing(r)
      return PluralRing{T}(r, basering)
   else
      n = Int(libSingular.rIsLPRing(r))
      if n > 0
         deg_bound = divexact(Int(libSingular.rVar(r)), n)
         return LPRing{T}(r, basering, deg_bound)
      else
         return PolyRing{T}(r, basering)
      end
   end
end

# for the translation of any return (tuple or non-tuple)
function convert_return(value::Vector, R = nothing)
    if !(length(value) > 0 && value[1] isa Bool)
        # value is a normal singular list
        return [convert_return(i, R) for i in value]
    elseif value[1]
        # value is a singular tuple: should only happen at the top level
        return [convert_return(value[i], R) for i in 2:length(value)]
    else
        # value is a normal singular value
        # singular lists here are behind pointers
        return convert_normal_value(value[2], value[3], R)
    end
end

function get_ring(arg_list)
    ring = nothing
    for i in arg_list
        current_ptr = nothing
        try
            current_ptr = i.ptr
        catch
            continue
        end
        if current_ptr isa poly
            return parent(i)
        elseif current_ptr isa ideal
            return parent(i).base_ring
        end
    end
    return ring
end

function prepare_argument(x::Vector{Int64})
    return Any[mapping_types_reversed[:INTVEC_CMD], libSingular.jl_array_to_intvec(x)], nothing
end

function prepare_argument(x::Matrix{Int64})
    return Any[mapping_types_reversed[:INTMAT_CMD], libSingular.jl_array_to_intmat(x)], nothing
end

function prepare_argument(x::Int64)
    return Any[mapping_types_reversed[:INT_CMD], Ptr{Cvoid}(x)], nothing
end

function prepare_argument(x::String)
    return Any[mapping_types_reversed[:STRING_CMD], libSingular.copy_string_to_void(x)], nothing
end

function prepare_argument(x::PolyRingUnion)
    GC.@preserve x new_ptr = libSingular.get_ring_ref(x.ptr)
    return Any[mapping_types_reversed[:RING_CMD], new_ptr], x
end

function prepare_argument(x::SPolyUnion)
    R = parent(x)
    GC.@preserve x R return (Any[mapping_types_reversed[:POLY_CMD],
                                 libSingular.copy_polyptr_to_void(x.ptr, R.ptr)],
                             R)
end

function prepare_argument(x::svector)
    R = parent(x).base_ring
    GC.@preserve x R return (Any[mapping_types_reversed[:VECTOR_CMD],
                                 libSingular.copy_polyptr_to_void(x.ptr, R.ptr)],
                             R)
end

function prepare_argument(x::sideal)
    R = parent(x).base_ring
    GC.@preserve x R return (Any[mapping_types_reversed[:IDEAL_CMD],
                                 libSingular.copy_idealptr_to_void(x.ptr, R.ptr)],
                             R)
end

function prepare_argument(x::sbigintmat)
    GC.@preserve x return (Any[mapping_types_reversed[:BIGINTMAT_CMD],
                                 libSingular.copy_bigintmatptr_to_void(x.ptr)],
                           false)
end

function prepare_argument(x::Union{
                          Matrix{ <: Union{Nemo.Integer, Nemo.ZZRingElem}},
                          Nemo.MatElem{ <: Union{Nemo.Integer, Nemo.ZZRingElem}},
                          Nemo.MatRingElem{ <: Union{Nemo.Integer, Nemo.ZZRingElem}}})
    return prepare_argument(sbigintmat(x))
end


function prepare_argument(x::Vector{Any}, R::PolyRingUnion)
    args = Vector{Any}()
    types = Vector{Int}()
    for i in x
       if typeof(i) == Vector{Any}
          p = prepare_argument(i, R)
       else
          p = prepare_argument(i)
       end
       push!(args, p[1][2])
       push!(types, p[1][1])
    end
    GC.@preserve x R return (Any[mapping_types_reversed[:LIST_CMD],
                                 libSingular.jl_array_to_void(args, types, R.ptr)],
                             R)
end

function prepare_argument(x::smodule)
    R = parent(x).base_ring
    GC.@preserve x R return (Any[mapping_types_reversed[:MODUL_CMD],
                                 libSingular.copy_idealptr_to_void(x.ptr, R.ptr)],
                             R)
end

function prepare_argument(x::sresolution)
    R = base_ring(x)
    GC.@preserve x R return (Any[mapping_types_reversed[:RESOLUTION_CMD],
                                 libSingular.create_syStrategy_data(x.ptr, R.ptr)],
                             R)
end

function prepare_argument(x::smatrix)
    R = base_ring(x)
    GC.@preserve x R return (Any[mapping_types_reversed[:MATRIX_CMD],
                                  libSingular.mp_Copy(x.ptr, R.ptr).cpp_object],
                              R)
end

function prepare_argument(x::BigInt)
    y = libSingular.number_ptr(x, libSingular.get_coeffs_BIGINT())
    return Any[mapping_types_reversed[:BIGINT_CMD], y.cpp_object], nothing
end

# errors
function prepare_argument(x::Vector)
   error("`intvec` may be passed in as Vector{Int}. All other vectors (`list` "*
          "in Singular) must be passed in as Vector{Any} along with an "*
          "explicit base ring in the first argument")
end

function prepare_argument(x::Any)
    :ptr in fieldnames(typeof(x)) || error("unrecognized argument $x")
    if x.ptr isa libSingular.number_ptr
        ptr = x.ptr
        ring_ptr = parent(x).ptr
        new_ptr = libSingular.n_Copy(ptr, ring_ptr)
        return Any[mapping_types_reversed[:NUMBER_CMD], new_ptr.cpp_object], nothing
    end
    error("unrecognized argument $x")
end

function low_level_caller_rng(lib::String, name::String, ring, args)
    libSingular.load_library(lib)
    arguments = Any[]
    for i in args
        if i isa Vector{Any}
           i, _ = prepare_argument(i, ring)
        else
           i, _ = prepare_argument(i)
        end
        push!(arguments, i)
    end
    return_value = libSingular.call_singular_library_procedure(name, ring.ptr, arguments)
    if libSingular.have_error()
      error(libSingular.get_and_clear_error())
    end
    return convert_return(return_value, ring)
end

function low_level_caller(lib::String, name::String, args)
    libSingular.load_library(lib)
    arguments = Any[]
    ring = nothing
    for i in args
        i, j = prepare_argument(i)
        push!(arguments, i)
        if j !== nothing
            ring = j
        end
    end
    ring_ptr = (ring === nothing) ? C_NULL : ring.ptr
    return_value = libSingular.call_singular_library_procedure(name, ring_ptr, arguments)
    if libSingular.have_error()
      error(libSingular.get_and_clear_error())
    end
    return convert_return(return_value, ring)
end

@doc raw"""
    lookup_library_symbol(package::String, name::String)

Attempt to look up a symbol in a particular Singular interpreter package and
return its value as a usable Singular.jl object. The package at the top level
is called "Top", and ring dependent objects are contained in their basering,
which is returned as a dictionary.

# Examples
```jldoctest
julia> Singular.call_interpreter("bigint a = 42;");

julia> a = Singular.lookup_library_symbol("Top", "a"); (a, typeof(a))
(42, BigInt)

julia> Singular.call_interpreter("ring r=0,(x,y,z),dp; poly f = (x+y)^2;");

julia> Singular.lookup_library_symbol("Top", "r")
2-element Vector{Any}:
 Singular polynomial ring (QQ),(x,y,z),(dp(3),C)
 Dict{Symbol, spoly{n_Q}}(:f => x^2 + 2*x*y + y^2)
```
"""
function lookup_library_symbol(pack::String, name::String)
    (err::Int, res) = libSingular.lookup_singular_library_symbol_wo_rng(pack, name)
    if err == 0
        return convert_return(res, nothing)
    elseif err == 1
        error("Singular symbol "*pack*"::"*name*" does not exist")
    else
        error("Singular package "*pack*" does not exist")
    end
end
