###############################################################################
#
#   SingularIntegerRing/n_Z
#
###############################################################################

type SingularIntegerRing <: Nemo.Ring
   ptr::libSingular.coeffs
   refcount::Int

   function SingularIntegerRing() 
      n_Z = @cxx n_Z
      d = new(libSingular.nInitChar(n_Z, Ptr{Void}(0)), 1)
      finalizer(d, _SingularIntegerRing_clear_fn)
      return d
   end
end

function _SingularIntegerRing_clear_fn(cf::SingularIntegerRing)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
end

type n_Z <: Nemo.RingElem
    ptr::libSingular.number

    function n_Z()
    	const c = SingularZZ.ptr
        z = new(libSingular.n_Init(0, c))
        parent(z).refcount += 1
        finalizer(z, _n_Z_clear_fn)
        return z
    end

    function n_Z(n::Int)
    	const c = SingularZZ.ptr
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
   _SingularIntegerRing_clear_fn(R)
   nothing
end

###############################################################################
#
#   SingularRationalField/n_Q
#
###############################################################################

type SingularRationalField <: Nemo.Field
   ptr::libSingular.coeffs
   refcount::Int

   function SingularRationalField() 
      n_Q = @cxx n_Q
      d = new(libSingular.nInitChar(n_Q, Ptr{Void}(0)), 1)
      finalizer(d, _SingularRationalField_clear_fn)
      return d
   end
end

function _SingularRationalField_clear_fn(cf::SingularRationalField)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

type n_Q <: Nemo.FieldElem
    ptr::libSingular.number

    function n_Q()
    	const c = SingularQQ.ptr
        z = new(libSingular.n_Init(0, c))
        parent(z).refcount += 1
        finalizer(z, _n_Q_clear_fn)
        return z
    end

    function n_Q(n::Int)
    	const c = SingularQQ.ptr
        z = new(libSingular.n_Init(n, c))
        parent(z).refcount += 1
        finalizer(z, _n_Q_clear_fn)
        return z
    end

    function n_Q(n::n_Z)
    	z = new(libSingular.nApplyMapFunc(n_Z_2_n_Q, n.ptr, SingularZZ.ptr, SingularQQ.ptr))
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
   _SingularRationalField_clear_fn(R)
   nothing
end

###############################################################################
#
#   SingularN_ZnRing/n_Zn
#
###############################################################################

type ZnmInfo
   n::BigInt
   exp::UInt
end

type SingularN_ZnRing <: Nemo.Ring
   ptr::libSingular.coeffs
   from_n_Z::Cxx.CppFptr
   to_n_Z::Cxx.CppFptr
   refcount::Int

   function SingularN_ZnRing(n::Int) 
      n_Zn = @cxx n_Zn
      ptr = libSingular.nInitChar(n_Zn, pointer_from_objref(ZnmInfo(BigInt(n), UInt(1))))
      d = new(ptr, libSingular.n_SetMap(SingularZZ.ptr, ptr), 
              libSingular.n_SetMap(ptr, SingularZZ.ptr), 1)
      finalizer(d, _SingularN_ZnRing_clear_fn)
      return d
   end
end

function _SingularN_ZnRing_clear_fn(cf::SingularN_ZnRing)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

type n_Zn <: Nemo.RingElem
    ptr::libSingular.number
    parent::SingularN_ZnRing

    function n_Zn(c::SingularN_ZnRing)
    	z = new(libSingular.n_Init(0, c.ptr))
        c.refcount += 1
        finalizer(z, _n_Zn_clear_fn)
        return z
    end

    function n_Zn(c::SingularN_ZnRing, n::Int)
    	z = new(libSingular.n_Init(n, c.ptr))
        c.refcount += 1
        finalizer(z, _n_Zn_clear_fn)
        return z
    end

    function n_Zn(c::SingularN_ZnRing, n::libSingular.number)
    	z = new(n)
        c.refcount += 1
        finalizer(z, _n_Zn_clear_fn)
        return z
    end
end

function _n_Zn_clear_fn(n::n_Zn)
   R = parent(n)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   _SingularN_ZnRing_clear_fn(R)
   nothing
end

###############################################################################
#
#   SingularN_ZpField/n_Zp
#
###############################################################################

type SingularN_ZpField <: Nemo.Field
   ptr::libSingular.coeffs
   from_n_Z::Cxx.CppFptr
   to_n_Z::Cxx.CppFptr
   refcount::Int

   function SingularN_ZpField(n::Int) 
      n_Zp = @cxx n_Zp
      ptr = libSingular.nInitChar(n_Zp, Ptr{Void}(n))
      d = new(ptr, libSingular.n_SetMap(SingularZZ.ptr, ptr), 
              libSingular.n_SetMap(ptr, SingularZZ.ptr), 1)
      finalizer(d, _SingularN_ZpField_clear_fn)
      return d
   end
end

function _SingularN_ZpField_clear_fn(cf::SingularN_ZpField)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

type n_Zp <: Nemo.FieldElem
    ptr::libSingular.number
    parent::SingularN_ZpField

    function n_Zp(c::SingularN_ZpField)
    	z = new(libSingular.n_Init(0, c.ptr))
        c.refcount += 1
        finalizer(z, _n_Zp_clear_fn)
        return z
    end

    function n_Zp(c::SingularN_ZpField, n::Int)
    	z = new(libSingular.n_Init(n, c.ptr))
        c.refcount += 1
        finalizer(z, _n_Zp_clear_fn)
        return z
    end

    function n_Zp(c::SingularN_ZpField, n::libSingular.number)
    	z = new(n)
        c.refcount += 1
        finalizer(z, _n_Zp_clear_fn)
        return z
    end
end

function _n_Zp_clear_fn(n::n_Zp)
   R = parent(n)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   _SingularN_ZpField_clear_fn(R)
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

type SingularN_GFField <: Nemo.Field
   ptr::libSingular.coeffs
   deg::Int
   from_n_Z::Cxx.CppFptr
   to_n_Z::Cxx.CppFptr
   refcount::Int

   function SingularN_GFField(p::Int, n::Int, S::Symbol) 
      n_GF = @cxx n_GF
      ptr = libSingular.nInitChar(n_GF, pointer_from_objref(GFInfo(Cint(p), Cint(n), pointer(Vector{UInt8}(string(S)*"\0")))))
      d = new(ptr, n, libSingular.n_SetMap(SingularZZ.ptr, ptr), 
              libSingular.n_SetMap(ptr, SingularZZ.ptr), 1)
      finalizer(d, _SingularN_GFField_clear_fn)
      return d
   end
end

function _SingularN_GFField_clear_fn(cf::SingularN_GFField)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
   nothing
end

type n_GF <: Nemo.FieldElem
    ptr::libSingular.number
    parent::SingularN_GFField

    function n_GF(c::SingularN_GFField)
    	z = new(libSingular.n_Init(0, c.ptr))
        c.refcount += 1
        finalizer(z, _n_GF_clear_fn)
        return z
    end

    function n_GF(c::SingularN_GFField, n::Int)
    	z = new(libSingular.n_Init(n, c.ptr))
        c.refcount += 1
        finalizer(z, _n_GF_clear_fn)
        return z
    end

    function n_GF(c::SingularN_GFField, n::libSingular.number)
    	z = new(n)
        c.refcount += 1
        finalizer(z, _n_GF_clear_fn)
        return z
    end
end

function _n_GF_clear_fn(n::n_GF)
   R = parent(n)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   _SingularN_GFField_clear_fn(R)
   nothing
end

###############################################################################
#
#   SingularCoefficientRing/n_unknown
#
###############################################################################

type CoefficientRing{T <: Nemo.RingElem} <: Nemo.Ring
   ptr::libSingular.coeffs
   base_ring::Nemo.Ring

   function CoefficientRing(R::Nemo.Ring)
      c = libSingular.register(R)
      ptr = pointer_from_objref(R)
      return new(libSingular.nInitChar(c, ptr), R)
   end
end

type n_unknown{T <: Nemo.RingElem} <: Nemo.RingElem
   ptr::libSingular.number
   parent::CoefficientRing{T}
end

