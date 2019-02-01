include("number/n_Z-test.jl")
include("number/n_Zp-test.jl")
include("number/n_Zn-test.jl")
include("number/n_GF-test.jl")
include("number/n_Q-test.jl")

function test_number()
   test_n_Z()
   test_n_Zp()
   test_n_Zn()
   test_n_GF()
   test_n_Q()
end
