###############################################################################
#
#   Integers/n_Z
#
###############################################################################

function get_n_Z()
   d = libSingular.nInitChar(libSingular.n_Z, Ptr{Void}(0))
end

const IntegersID = Dict{Symbol, Ring}()

type Integers <: Ring
   ptr::libSingular.coeffs
   refcount::Int

   function Integers()
      if haskey(IntegersID, :ZZ)
         d = IntegersID[:ZZ]::Integers
      else
         d = new()
         IntegersID[:ZZ] = d
         finalizer(d, _Integers_clear_fn)
      end
      return d
   end
end

function _Integers_clear_fn(cf::Integers)
   cf.refcount -= 1
   if cf.refcount == 0
      libSingular.nKillChar(cf.ptr)
   end
end

