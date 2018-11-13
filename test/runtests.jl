using Singular

import AbstractAlgebra
import Nemo

if VERSION < v"0.7.0-DEV.2004"
   using Base.Test
else
   using Test
end

include("../test/number-test.jl")
include("../test/poly-test.jl")
include("../test/ideal-test.jl")
include("../test/matrix-test.jl")
include("../test/resolution-test.jl")
include("../test/module-test.jl")
# include("../test/libsingular-test.jl")

test_number()
test_poly()
test_ideal()
test_matrix()
test_resolution()
test_module()
# test_libsingular()

