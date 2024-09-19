abstract type Set <: AbstractAlgebra.Set end

abstract type Ring <: AbstractAlgebra.Ring end

abstract type Field <: AbstractAlgebra.Field end

abstract type Module{T <: AbstractAlgebra.RingElem} <: AbstractAlgebra.Module{T} end

abstract type Map <: AbstractAlgebra.SetMap end

abstract type AbstractAlgebraHomomorphism <: Map end

