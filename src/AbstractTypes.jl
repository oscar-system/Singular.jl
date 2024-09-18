abstract type Set <: AbstractAlgebra.Set end

abstract type NCRing <: AbstractAlgebra.NCRing end

abstract type Ring <: AbstractAlgebra.Ring end

abstract type Field <: AbstractAlgebra.Field end

abstract type Module{T <: AbstractAlgebra.NCRingElem} <: AbstractAlgebra.Module{T} end

abstract type Map <: AbstractAlgebra.SetMap end

abstract type AbstractAlgebraHomomorphism <: Map end

