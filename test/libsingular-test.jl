include("libsingular/nemo-test.jl")
include("libsingular/Puiseuxexpansions-test.jl")

@testset "libsingular.number_layout" begin
   # Singular.jl reads n_Z and n_Q coefficients directly out of Singular's
   # memory; this is the check that __init__ already ran
   @test Singular.libSingular.check_number_layout(ZZ.ptr, QQ.ptr) === nothing
end
