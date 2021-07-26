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

# Finally the really fun part: When coverting everything back to something for
# julia, the following both produce the same julia result [a, b]:
#     - (a,b): a tuple of length 2 returned by a procedure
#     - list(a,b): a list of length 2 returned by a procedure
# Therefore, the user has to know how to interpret the result.

casting_functions = nothing

function create_casting_functions()
    return Dict(
        mapping_types_reversed[:NUMBER_CMD] =>
            function (vptr, R)
                cast = libSingular.NUMBER_CMD_CASTER(vptr)
                # TODO hmm, the number should not be put into the poly ring
                return R(cast)
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
                return cast
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
                return cast
            end
        ,
        mapping_types_reversed[:INTMAT_CMD] =>
            function (vptr, R)
                cast = libSingular.INTMAT_CMD_CASTER(vptr)
                return cast
            end
        ,
        mapping_types_reversed[:BIGINT_CMD] =>
            function (vptr, R)
                error("unable to return a bigint from a singular library procedure")
#                cast = libSingular.NUMBER_CMD_CASTER(vptr)
#                return libSingular.n_GetMPZ(cast, libSingular.get_coeffs_BIGINT())
            end
        ,
        mapping_types_reversed[:BIGINTMAT_CMD] =>
            function (vptr, R)
                error("unable to return a bigintmat from a singular library procedure")
                cast = libSingular.BIGINTMAT_CMD_CASTER(vptr)
                # TODO this leak cannot possibly be useful
                return cast
            end
        ,
        mapping_types_reversed[:MAP_CMD] =>
            function (vptr, R)
                @warn "returning a map from a singular library procedure as an ideal" maxlog=1
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
                    error("unrecognized return with singular type number $typ")
                end
    return mapper(valueptr, R)
end

function convert_ring_content(value_list, R)
    return Dict(i[2] => convert_normal_value(i[3], i[1], R) for i in value_list)
end

# take ownership of r
function create_ring_from_singular_ring(r::libSingular.ring_ptr)
   c = libSingular.rCoeffPtr(r)
   if libSingular.nCoeff_is_Q(c)
      return PolyRing{n_Q}(r, QQ)
   elseif libSingular.nCoeff_is_Zp(c)
      p = Int(libSingular.n_GetChar(c))
      return PolyRing{n_Zp}(r, N_ZpField(p))
   elseif libSingular.nCoeff_is_GF(c)
      p = Int(libSingular.n_GetChar(c))
      q = Int(libSingular.nfCharQ(c))
      d = round(Int, log(p, q))
      s = Symbol(libSingular.n_ParameterName(0, c))
      return PolyRing{n_GF}(r, N_GField(p, d, s))
   elseif libSingular.nCoeff_is_transExt(c)
      p = Int(libSingular.n_GetChar(c))
      npars = libSingular.n_NumberOfParameters(c)
      S = [Symbol(libSingular.n_ParameterName(i-1, c)) for i in 1:npars]
      F = N_FField(iszero(p) ? QQ : N_ZpField(p), S)
      return PolyRing{n_transExt}(r, F)
   elseif libSingular.nCoeff_is_algExt(c)
      # first create the univariate transcendental extension
      p = Int(libSingular.n_GetChar(c))
      @assert libSingular.n_NumberOfParameters(c) == 1
      S = [Symbol(libSingular.n_ParameterName(0, c))]
      F = N_FField(iszero(p) ? QQ : N_ZpField(p), S)
      # now create the extension
      minpoly = F(libSingular.algExt_GetMinpoly(c, F.ptr))
      K = N_AlgExtField(libSingular.nCopyCoeff(c), minpoly)
      return PolyRing{n_algExt}(r, K)
   else
      basering = N_UnknownSingularCoefficientRing(libSingular.nCopyCoeff(c))
      return PolyRing{n_unknownsingularcoefficient}(r, basering)
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

function prepare_argument(x::Array{Int64, 1})
    return Any[mapping_types_reversed[:INTVEC_CMD], libSingular.jl_array_to_intvec(x)], nothing
end

function prepare_argument(x::Array{Int64, 2})
    return Any[mapping_types_reversed[:INTMAT_CMD], libSingular.jl_array_to_intmat(x)], nothing
end

function prepare_argument(x::Int64)
    return Any[mapping_types_reversed[:INT_CMD], Ptr{Cvoid}(x)], nothing
end

function prepare_argument(x::String)
    return Any[mapping_types_reversed[:STRING_CMD], libSingular.copy_string_to_void(x)], nothing
end

function prepare_argument(x::PolyRing)
    GC.@preserve x new_ptr = libSingular.get_ring_ref(x.ptr)
    return Any[mapping_types_reversed[:RING_CMD], new_ptr], x
end

function prepare_argument(x::spoly)
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

function prepare_argument(x::Array{Any, 1}, R::PolyRing)
    args = Array{Any, 1}()
    types = Array{Int, 1}()
    for i in x
       if typeof(i) == Array{Any, 1}
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

function prepare_argument(x::Any)
    :ptr in fieldnames(typeof(x)) || error("unrecognized argument $x")
    if x.ptr isa libSingular.number_ptr
        ptr = x.ptr
        rng = parent(x)
        new_ptr = libSingular.n_Copy(ptr, rng.ptr)
        return Any[mapping_types_reversed[:NUMBER_CMD], new_ptr.cpp_object], nothing
    elseif x.ptr isa libSingular.__mpz_struct
        return Any[mapping_types_reversed[:BIGINT_CMD], x.ptr.cpp_object], nothing
    elseif x.ptr isa libSingular.map_ptr
        return Any[mapping_types_reversed[:MAP_CMD], x.ptr.cpp_object], nothing
    elseif x.ptr isa libSingular.bigintmat
        return Any[mapping_types_reversed[:BIGINTMAT_CMD], x.ptr.cpp_object], nothing
    end
    error("unrecognized argument $x")
end

function low_level_caller_rng(lib::String, name::String, ring, args...)
    libSingular.load_library(lib)
    arguments = Array{Any, 1}()
    for i in args
       if typeof(i) == Array{Any, 1} 
          push!(arguments, prepare_argument(i, ring))
       else
          push!(arguments, prepare_argument(i))
       end
    end
    arguments = Any[i for (i, j) in arguments]
    return_value = libSingular.call_singular_library_procedure(name, ring.ptr, arguments)
    return convert_return(return_value, ring)
end

function low_level_caller(lib::String, name::String, args...)
    libSingular.load_library(lib)
    arguments = [prepare_argument(i) for i in args]
    rng = nothing
    for (i, j) in arguments
        if j != nothing
            rng = j
        end
    end
    arguments = Any[i for (i, j) in arguments]
    return_values = nothing
    rng_ptr = (rng == nothing) ? C_NULL : rng.ptr
    return_value = libSingular.call_singular_library_procedure(name, rng_ptr, arguments)
    return convert_return(return_value, rng)
end
