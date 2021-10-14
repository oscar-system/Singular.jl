using Singular

using Random: Random, MersenneTwister
const rng = MersenneTwister()

using RandomExtensions: make

import AbstractAlgebra
import Nemo

using Test

include("../test/number-test.jl")
include("../test/poly-test.jl")
include("../test/exterior-test.jl")
include("../test/weyl-test.jl")
include("../test/ideal-test.jl")
include("../test/map-test.jl")
include("../test/matrix-test.jl")
include("../test/resolution-test.jl")
include("../test/module-test.jl")
include("../test/call_interpreter-test.jl")
include("../test/caller-test.jl")
include("../test/libsingular-test.jl")

