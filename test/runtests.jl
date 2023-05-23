using Singular

using Singular.Random: Random, MersenneTwister
const rng = MersenneTwister()

using Singular.RandomExtensions: make

import Singular.AbstractAlgebra
import Singular.Nemo

using Test

include("number-test.jl")
include("poly-test.jl")
include("ideal-test.jl")
include("map-test.jl")
include("matrix-test.jl")
include("resolution-test.jl")
include("module-test.jl")
include("call_interpreter-test.jl")
include("caller-test.jl")
include("libsingular-test.jl")
