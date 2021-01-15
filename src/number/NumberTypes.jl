###############################################################################
#
#   Finalizers
#
###############################################################################

function _Ring_finalizer(cf)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

function _Elem_finalizer(n)
   R = parent(n)
   libSingular.n_Delete(n.ptr, R.ptr)
   _Ring_finalizer(R)
   nothing
end

###############################################################################
#
#   Integers/n_Z
#
###############################################################################

using CxxWrap

const IntegersID = Dict{Symbol, Ring}()

mutable struct Integers <: Ring
   ptr::libSingular.coeffs_ptr
   refcount::Int

   function Integers()
      if haskey(IntegersID, :ZZ)
         d = IntegersID[:ZZ]::Integers
      else
         ptr = libSingular.nInitChar(libSingular.n_Z, Ptr{Nothing}(0))
         d = new(ptr, 0)
         IntegersID[:ZZ] = d
         finalizer(_Ring_finalizer, d)
      end
      return d
   end
end

mutable struct n_Z <: Nemo.RingElem
    ptr::libSingular.number_ptr

    function n_Z(n::libSingular.number_ptr)
        z = new(n)
        parent(z).refcount += 1
        finalizer(_Elem_finalizer, z)
        return z
    end
end

n_Z(n::Int = 0) = n_Z(libSingular.n_Init(n, ZZ.ptr))


###############################################################################
#
#   Rationals/n_Q
#
###############################################################################

const RationalsID = Dict{Symbol, Field}()

mutable struct Rationals <: Field
   ptr::libSingular.coeffs_ptr
   refcount::Int

   function Rationals()
      if haskey(RationalsID, :QQ)
         d = RationalsID[:QQ]::Rationals
      else
         ptr = libSingular.nInitChar(libSingular.n_Q, Ptr{Nothing}(0))
         d = new(ptr, 0)
         RationalsID[:QQ] = d
         finalizer(_Ring_finalizer, d)
      end
      return d
   end
end

mutable struct n_Q <: Nemo.FieldElem
    ptr::libSingular.number_ptr

    function n_Q(n::libSingular.number_ptr)
        z = new(n)
        parent(z).refcount += 1
        finalizer(_Elem_finalizer, z)
        return z
    end
end

n_Q(n::Int = 0) = n_Q(libSingular.n_Init(n, QQ.ptr))
n_Q(n::n_Z) = n_Q(libSingular.nApplyMapFunc(n_Z_2_n_Q, n.ptr, ZZ.ptr, QQ.ptr))


###############################################################################
#
#   N_ZnRing/n_Zn
#
###############################################################################

mutable struct ZnmInfo
   n::BigInt
   exp::UInt
end

const N_ZnRingID = Dict{Int, Ring}()

mutable struct N_ZnRing <: Ring
   ptr::libSingular.coeffs_ptr
   from_n_Z::Ptr{Nothing}
   to_n_Z::Ptr{Nothing}
   refcount::Int

   function N_ZnRing(n::Int)
      if haskey(N_ZnRingID, n)
         d = N_ZnRingID[n]::N_ZnRing
      else
         ptr = libSingular.nInitChar(libSingular.n_Zn, pointer_from_objref(ZnmInfo(BigInt(n), UInt(1))))
         d = new(ptr, libSingular.n_SetMap(ZZ.ptr, ptr),
              libSingular.n_SetMap(ptr, ZZ.ptr), 1)
         N_ZnRingID[n] = d
         finalizer(_Ring_finalizer, d)
      end
      return d
   end
end


mutable struct n_Zn <: Nemo.RingElem
    ptr::libSingular.number_ptr
    parent::N_ZnRing

    function n_Zn(c::N_ZnRing, n::libSingular.number_ptr)
        z = new(n, c)
        c.refcount += 1
        finalizer(_Elem_finalizer, z)
        return z
    end
end

n_Zn(c::N_ZnRing, n::Int = 0) = n_Zn(c, libSingular.n_Init(n, c.ptr))



###############################################################################
#
#   N_ZpField/n_Zp
#
###############################################################################

const N_ZpFieldID = Dict{Int, Field}()

mutable struct N_ZpField <: Field
   ptr::libSingular.coeffs_ptr
   from_n_Z::Ptr{Nothing}
   to_n_Z::Ptr{Nothing}
   refcount::Int

   function N_ZpField(n::Int)
      if haskey(N_ZpFieldID, n)
         d = N_ZpFieldID[n]::N_ZpField
      else
         ptr = libSingular.nInitChar(libSingular.n_Zp, Ptr{Nothing}(n))
         d = new(ptr, libSingular.n_SetMap(ZZ.ptr, ptr),
              libSingular.n_SetMap(ptr, ZZ.ptr), 1)
         N_ZpFieldID[n] = d
         finalizer(_Ring_finalizer, d)
      end
      return d
   end
end

mutable struct n_Zp <: Nemo.FieldElem
    ptr::libSingular.number_ptr
    parent::N_ZpField

    function n_Zp(c::N_ZpField, n::libSingular.number_ptr)
        z = new(n, c)
        c.refcount += 1
        finalizer(_Elem_finalizer, z)
        return z
    end
end

n_Zp(c::N_ZpField, n::Int = 0) = n_Zp(c, libSingular.n_Init(n, c.ptr))


###############################################################################
#
#   SingularFiniteField/n_GF
#
###############################################################################

mutable struct GFInfo
   p::Cint
   n::Cint
   s::Ptr{UInt8}
end

const N_GFieldID = Dict{Tuple{Int, Int, Symbol}, Field}()

mutable struct N_GField <: Field
   ptr::libSingular.coeffs_ptr
   deg::Int
   from_n_Z::Ptr{Nothing}
   to_n_Z::Ptr{Nothing}
   refcount::Int
   S::Symbol

   function N_GField(p::Int, n::Int, S::Symbol)
      if haskey(N_GFieldID, (p, n, S))
         d = N_GFieldID[p, n, S]::N_GField
      else
         gfinfo = GFInfo(Cint(p), Cint(n), pointer(Base.Vector{UInt8}(string(S)*"\0")))
         GC.@preserve gfinfo begin
         ptr = libSingular.nInitChar(libSingular.n_GF, pointer_from_objref(gfinfo))
         d = new(ptr, n, libSingular.n_SetMap(ZZ.ptr, ptr),
              libSingular.n_SetMap(ptr, ZZ.ptr), 1, S)
         end
         N_GFieldID[p, n, S] = d
         finalizer(_Ring_finalizer, d)
      end
      return d
   end
end

mutable struct n_GF <: Nemo.FieldElem
    ptr::libSingular.number_ptr
    parent::N_GField

    function n_GF(c::N_GField, n::libSingular.number_ptr)
        z = new(n, c)
        c.refcount += 1
        finalizer(_Elem_finalizer, z)
        return z
    end
end

n_GF(c::N_GField, n::Int = 0) = n_GF(c, libSingular.n_Init(n, c.ptr))


###############################################################################
#
#   SingularTranscendentalExtensionField/n_transExt
#
###############################################################################

const N_FFieldID = Dict{Tuple{Field, Array{Symbol, 1}}, Field}()

mutable struct N_FField <: Field
   ptr::libSingular.coeffs_ptr
   base_ring::Field
   refcount::Int

   function N_FField(F::Field, S::Array{Symbol, 1})
      if haskey(N_FFieldID, (F, S))
         d = N_FFieldID[F, S]::N_FField
      else
         v = [pointer(Base.Vector{UInt8}(string(str)*"\0")) for str in S]
         cf = libSingular.nCopyCoeff(F.ptr)
         ptr = libSingular.transExt_helper(cf, v)
         d = new(ptr, F, 1)
         N_FFieldID[F, S] = d
         finalizer(_Ring_finalizer, d)
      end
      return d
   end
end

mutable struct n_transExt <: Nemo.FieldElem
    ptr::libSingular.number_ptr
    parent::N_FField

    function n_transExt(c::N_FField, n::libSingular.number_ptr)
        z = new(n, c)
        c.refcount += 1
        finalizer(_Elem_finalizer, z)
        return z
    end
end

n_transExt(c::N_FField, n::Int = 0) = n_transExt(c, libSingular.n_Init(n, c.ptr))


###############################################################################
#
#   SingularCoefficientRing/n_unknown
#
###############################################################################

const CoeffRingID = Dict{Nemo.Ring, Ring}()

mutable struct CoefficientRing{T <: Nemo.RingElem} <: Ring
   ptr::libSingular.coeffs_ptr
   base_ring::Nemo.Ring

   function CoefficientRing{T}(R::Nemo.Ring, cached::Bool=true) where T
      if haskey(CoeffRingID, R)
         return CoeffRingID[R]::CoefficientRing{T}
      else
         c = libSingular.register(R)
         ptr = pointer_from_objref(R)
         z = new(libSingular.nInitChar(c, ptr), R)
         if cached
           CoeffRingID[R] = z
         end
         return z
      end
   end
end

mutable struct n_unknown{T <: Nemo.RingElem} <: Nemo.RingElem
   ptr::libSingular.number_ptr
   parent::CoefficientRing{T}

end

###############################################################################
#
#   N_UnknownSingularCoefficientRing
#   singular coefficient rings with no identifiable julia object
#
###############################################################################

mutable struct N_UnknownSingularCoefficientRing <: Ring
   ptr::libSingular.coeffs_ptr
   refcount::Int

   # take ownership of the pointer
   function N_UnknownSingularCoefficientRing(ptr::libSingular.coeffs_ptr)
      a = new(ptr, 1)
      finalizer(_Ring_finalizer, a)
      return a
   end
end

mutable struct n_unknownsingularcoefficient <: Nemo.RingElem
   ptr::libSingular.number_ptr
   parent::N_UnknownSingularCoefficientRing

   function n_unknownsingularcoefficient(R::N_UnknownSingularCoefficientRing, n::libSingular.number_ptr)
      a = new(n, R)
      finalizer(_Elem_finalizer, a)
      return a
   end
end
