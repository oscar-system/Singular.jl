###############################################################################
#
#   Integers/n_Z
#
###############################################################################

function get_n_Z()
   n_Z = @cxx n_Z
   d = libSingular.nInitChar(n_Z, Ptr{Void}(0))
end

type Integers <: Ring
   ptr::libSingular.coeffs
   refcount::Int

   function Integers() 
      d = new()
      finalizer(d, _Integers_clear_fn)
      return d
   end
end

function _Integers_clear_fn(cf::Integers)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
end

type n_Z <: Nemo.RingElem
    ptr::libSingular.number

    function n_Z()
    	const c = ZZ.ptr
        z = new(libSingular.n_Init(0, c))
        parent(z).refcount += 1
        finalizer(z, _n_Z_clear_fn)
        return z
    end

    function n_Z(n::Int)
    	const c = ZZ.ptr
        z = new(libSingular.n_Init(n, c))
        parent(z).refcount += 1
        finalizer(z, _n_Z_clear_fn)
        return z
    end

    function n_Z(n::libSingular.number)
    	z = new(n)
        parent(z).refcount += 1
        finalizer(z, _n_Z_clear_fn)
        return z
    end
end

function _n_Z_clear_fn(n::n_Z)
   R = parent(n)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   _Integers_clear_fn(R)
   nothing
end

###############################################################################
#
#   RationalField/n_Q
#
###############################################################################

function get_n_Q()
   n_Q = @cxx n_Q
   d = libSingular.nInitChar(n_Q, Ptr{Void}(0))
end

type RationalField <: Field
   ptr::libSingular.coeffs
   refcount::Int

   function RationalField() 
      d = new()
      finalizer(d, _RationalField_clear_fn)
      return d
   end
end

function _RationalField_clear_fn(cf::RationalField)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

type n_Q <: Nemo.FieldElem
    ptr::libSingular.number

    function n_Q()
    	const c = QQ.ptr
        z = new(libSingular.n_Init(0, c))
        parent(z).refcount += 1
        finalizer(z, _n_Q_clear_fn)
        return z
    end

    function n_Q(n::Int)
    	const c = QQ.ptr
        z = new(libSingular.n_Init(n, c))
        parent(z).refcount += 1
        finalizer(z, _n_Q_clear_fn)
        return z
    end

    function n_Q(n::n_Z)
    	z = new(libSingular.nApplyMapFunc(n_Z_2_n_Q, n.ptr, ZZ.ptr, QQ.ptr))
        parent(z).refcount += 1
        finalizer(z, _n_Q_clear_fn)
        return z
    end

    function n_Q(n::libSingular.number)
    	z = new(n)
        parent(z).refcount += 1
        finalizer(z, _n_Q_clear_fn)
        return z
    end
end

function _n_Q_clear_fn(n::n_Q)
   R = parent(n)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   _RationalField_clear_fn(R)
   nothing
end

###############################################################################
#
#   N_ZnRing/n_Zn
#
###############################################################################

type ZnmInfo
   n::BigInt
   exp::UInt
end

type N_ZnRing <: Ring
   ptr::libSingular.coeffs
   from_n_Z::Cxx.CppFptr
   to_n_Z::Cxx.CppFptr
   refcount::Int

   function N_ZnRing(n::Int) 
      n_Zn = @cxx n_Zn
      ptr = libSingular.nInitChar(n_Zn, pointer_from_objref(ZnmInfo(BigInt(n), UInt(1))))
      d = new(ptr, libSingular.n_SetMap(ZZ.ptr, ptr), 
              libSingular.n_SetMap(ptr, ZZ.ptr), 1)
      finalizer(d, _N_ZnRing_clear_fn)
      return d
   end
end

function _N_ZnRing_clear_fn(cf::N_ZnRing)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

type n_Zn <: Nemo.RingElem
    ptr::libSingular.number
    parent::N_ZnRing

    function n_Zn(c::N_ZnRing)
    	z = new(libSingular.n_Init(0, c.ptr))
        c.refcount += 1
        finalizer(z, _n_Zn_clear_fn)
        return z
    end

    function n_Zn(c::N_ZnRing, n::Int)
    	z = new(libSingular.n_Init(n, c.ptr))
        c.refcount += 1
        finalizer(z, _n_Zn_clear_fn)
        return z
    end

    function n_Zn(c::N_ZnRing, n::libSingular.number)
    	z = new(n)
        c.refcount += 1
        finalizer(z, _n_Zn_clear_fn)
        return z
    end
end

function _n_Zn_clear_fn(n::n_Zn)
   R = parent(n)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   _N_ZnRing_clear_fn(R)
   nothing
end

###############################################################################
#
#   N_ZpField/n_Zp
#
###############################################################################

type N_ZpField <: Field
   ptr::libSingular.coeffs
   from_n_Z::Cxx.CppFptr
   to_n_Z::Cxx.CppFptr
   refcount::Int

   function N_ZpField(n::Int) 
      n_Zp = @cxx n_Zp
      ptr = libSingular.nInitChar(n_Zp, Ptr{Void}(n))
      d = new(ptr, libSingular.n_SetMap(ZZ.ptr, ptr), 
              libSingular.n_SetMap(ptr, ZZ.ptr), 1)
      finalizer(d, _N_ZpField_clear_fn)
      return d
   end
end

function _N_ZpField_clear_fn(cf::N_ZpField)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

type n_Zp <: Nemo.FieldElem
    ptr::libSingular.number
    parent::N_ZpField

    function n_Zp(c::N_ZpField)
    	z = new(libSingular.n_Init(0, c.ptr))
        c.refcount += 1
        finalizer(z, _n_Zp_clear_fn)
        return z
    end

    function n_Zp(c::N_ZpField, n::Int)
    	z = new(libSingular.n_Init(n, c.ptr))
        c.refcount += 1
        finalizer(z, _n_Zp_clear_fn)
        return z
    end

    function n_Zp(c::N_ZpField, n::libSingular.number)
    	z = new(n)
        c.refcount += 1
        finalizer(z, _n_Zp_clear_fn)
        return z
    end
end

function _n_Zp_clear_fn(n::n_Zp)
   R = parent(n)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   _N_ZpField_clear_fn(R)
   nothing
end

###############################################################################
#
#   SingularFiniteField/n_GF
#
###############################################################################

type GFInfo
   p::Cint
   n::Cint
   s::Ptr{UInt8}
end

type N_GField <: Field
   ptr::libSingular.coeffs
   deg::Int
   from_n_Z::Cxx.CppFptr
   to_n_Z::Cxx.CppFptr
   refcount::Int

   function N_GField(p::Int, n::Int, S::Symbol) 
      n_GF = @cxx n_GF
      ptr = libSingular.nInitChar(n_GF, pointer_from_objref(GFInfo(Cint(p), Cint(n), pointer(Base.Vector{UInt8}(string(S)*"\0")))))
      d = new(ptr, n, libSingular.n_SetMap(ZZ.ptr, ptr), 
              libSingular.n_SetMap(ptr, ZZ.ptr), 1)
      finalizer(d, _N_GField_clear_fn)
      return d
   end
end

function _N_GField_clear_fn(cf::N_GField)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

type n_GF <: Nemo.FieldElem
    ptr::libSingular.number
    parent::N_GField

    function n_GF(c::N_GField)
    	z = new(libSingular.n_Init(0, c.ptr))
        c.refcount += 1
        finalizer(z, _n_GF_clear_fn)
        return z
    end

    function n_GF(c::N_GField, n::Int)
    	z = new(libSingular.n_Init(n, c.ptr))
        c.refcount += 1
        finalizer(z, _n_GF_clear_fn)
        return z
    end

    function n_GF(c::N_GField, n::libSingular.number)
    	z = new(n)
        c.refcount += 1
        finalizer(z, _n_GF_clear_fn)
        return z
    end
end

function _n_GF_clear_fn(n::n_GF)
   R = parent(n)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   _N_GField_clear_fn(R)
   nothing
end

###############################################################################
#
#   SingularCoefficientRing/n_unknown
#
###############################################################################

const CoeffRingID = Dict{Nemo.Ring, Nemo.Ring}()

type CoefficientRing{T <: Nemo.RingElem} <: Ring
   ptr::libSingular.coeffs
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

type n_unknown{T <: Nemo.RingElem} <: Nemo.RingElem
   ptr::libSingular.number
   parent::CoefficientRing{T}
end

