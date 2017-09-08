abstract type Set <: Nemo.Set end

abstract type Ring <: Nemo.Ring end

abstract type Field <: Nemo.Field end

abstract type Module{T <: Nemo.RingElem} <: Nemo.Module{T} end