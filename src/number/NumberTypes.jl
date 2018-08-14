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
#   Rationals/n_Q
#
###############################################################################

function get_n_Q()
   d = libSingular.nInitChar(libSingular.n_Q, Ptr{Void}(0))
end

const RationalsID = Dict{Symbol, Field}()

type Rationals <: Field
   ptr::libSingular.coeffs
   refcount::Int

   function Rationals()
      if haskey(RationalsID, :QQ)
         d = RationalsID[:QQ]::Rationals
      else
         d = new()
         RationalsID[:QQ] = d
         finalizer(d, _Rationals_clear_fn)
      end
      return d
   end
end

function _Rationals_clear_fn(cf::Rationals)
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
   _Rationals_clear_fn(R)
   nothing
end

