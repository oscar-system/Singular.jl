export CoeffsField

###############################################################################
#
#   Singular coefficient fields
#
###############################################################################

const CoeffsFieldID = ObjectIdDict() # All CoeffsFields are to be unique
const CoeffsPtrID = ObjectIdDict()

type CoeffsField <: Nemo.Field
   ptr::libSingular.coeffs

   function CoeffsField(nt::libSingular.n_coeffType) 
      par = Ptr{Void}(0)
      try
         return CoeffsFieldID[nt, par]
      catch
         ptr = libSingular.nInitChar(nt, par)
         #(ptr == libSingular.coeffs(0)) && error("Singular coeffs.domain construction failure")
         try
            cf = CoeffsPtrID[ptr]
	        libSingular.nKillChar(ptr)
	        return cf
         catch
            d = CoeffsFieldID[nt, par] = CoeffsPtrID[ptr] = new(ptr)
            finalizer(d, _CoeffsField_clear_fn)
            return d
         end
      end
   end
end

function _CoeffsField_clear_fn(cf::CoeffsField)
   libSingular.nKillChar(cf.ptr)
end

###############################################################################
#
#   SingularIntegerRing/n_Z
#
###############################################################################

type SingularIntegerRing <: Nemo.Ring
   ptr::libSingular.coeffs

   function SingularIntegerRing() 
      n_Z = @cxx n_Z
      d = new(libSingular.nInitChar(n_Z, Ptr{Void}(0)))
      finalizer(d, _SingularIntegerRing_clear_fn)
      return d
   end
end

function _SingularIntegerRing_clear_fn(cf::SingularIntegerRing)
   libSingular.nKillChar(cf.ptr)
end

type n_Z <: Nemo.RingElem
    ptr::libSingular.number

    function n_Z()
    	const c = SingularZZ.ptr
        z = new(libSingular.n_Init(0, c))
        finalizer(z, _n_Z_clear_fn)
        return z
    end

    function n_Z(n::Int)
    	const c = SingularZZ.ptr
        z = new(libSingular.n_Init(n, c))
        finalizer(z, _n_Z_clear_fn)
        return z
    end

    function n_Z(n::libSingular.number)
    	z = new(n)
        finalizer(z, _n_Z_clear_fn)
        return z
    end
end

function _n_Z_clear_fn(n::n_Z)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
end

###############################################################################
#
#   SingularRationalField/n_Q
#
###############################################################################

type SingularRationalField <: Nemo.Field
   ptr::libSingular.coeffs

   function SingularRationalField() 
      n_Q = @cxx n_Q
      d = new(libSingular.nInitChar(n_Q, Ptr{Void}(0)))
      finalizer(d, _SingularRationalField_clear_fn)
      return d
   end
end

function _SingularRationalField_clear_fn(cf::SingularRationalField)
   libSingular.nKillChar(cf.ptr)
end

type n_Q <: Nemo.FieldElem
    ptr::libSingular.number

    function n_Q()
    	const c = SingularQQ.ptr
        z = new(libSingular.n_Init(0, c))
        finalizer(z, _n_Q_clear_fn)
        return z
    end

    function n_Q(n::Int)
    	const c = SingularQQ.ptr
        z = new(libSingular.n_Init(n, c))
        finalizer(z, _n_Q_clear_fn)
        return z
    end

    function n_Q(n::n_Z)
    	z = new(libSingular.nApplyMapFunc(n_Z_2_n_Q, p, SingularZZ.ptr, SingularQQ.ptr))
        finalizer(z, _n_Q_clear_fn)
        return z
    end

    function n_Q(n::libSingular.number)
    	z = new(n)
        finalizer(z, _n_Q_clear_fn)
        return z
    end
end

function _n_Q_clear_fn(n::n_Q)
   libSingular.n_Delete(n.ptr, parent(n).ptr)
   nothing
end
