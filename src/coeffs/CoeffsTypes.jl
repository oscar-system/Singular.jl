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
#   SingularRationalField/SingularQQElem
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

type SingularQQElem <: Nemo.FieldElem
    ptr::libSingular.number

    function SingularQQElem()
    	const c = SingularQQ.ptr
        z = new(libSingular.n_Init(0, c))
        finalizer(z, _SingularQQElem_clear_fn)
        return z
    end

    function SingularQQElem(n::Int)
    	const c = SingularQQ.ptr
        z = new(libSingular.n_Init(n, c))
        finalizer(z, _SingularQQElem_clear_fn)
        return z
    end

    function SingularQQElem(n::libSingular.number)
    	z = new(n)
        finalizer(z, _SingularQQElem_clear_fn)
        return z
    end
end

function _SingularQQElem_clear_fn(n::SingularQQElem)
   c = parent(n)
   cf = c.ptr
   p = n.ptr

   # if libSingular.nCoeff_has_simple_Alloc(cf) || (p == number(0))
   #    nothing
   # end

   libSingular.n_Delete(p, cf)
end
