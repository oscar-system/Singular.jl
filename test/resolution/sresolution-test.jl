function test_sresolution_constructors()
   print("sresolution.constructors...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x*y + 1, x^2 + 1)
   F = fres(std(I), 3)
   S = parent(F)

   @test elem_type(S) == sresolution{spoly{n_Q}}
   @test elem_type(ResolutionSet{spoly{n_Q}}) == sresolution{spoly{n_Q}}
   @test parent_type(sresolution{spoly{n_Q}}) == ResolutionSet{spoly{n_Q}}

   @test base_ring(S) == R
   @test base_ring(F) == R

   @test typeof(S) <: AbstractAlgebra.Set
   @test typeof(F) <: AbstractAlgebra.SetElem

   @test isa(F, sresolution)

   println("PASS")
end

function test_sresolution_manipulation()
   print("sresolution.manipulation...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x*y + 1, x^2 + 1)
   F = fres(std(I), 1)

   @test length(F) == 1

   F = fres(std(I), 0)

   @test length(F) == 2

   F = fres(std(I), 2)

   @test length(F) == 2

   F = fres(std(I), 3)

   @test length(F) == 2

   @test isa(F[1], smodule)
   @test isa(F[2], smodule)
   @test isa(F[3], smodule)

   G = deepcopy(F)

   @test isa(G, sresolution)
   @test parent(F) === parent(G)

   println("PASS")
end

#= function test_sresolution_betti()
   print("sresolution.betti...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x*y + 1, x^2 + 1)
   F = fres(std(I), 3)
   M = minres(F)

   B = betti(M)

   @test size(B) == (2, 3)

   @test B[1, 1] == 1 && B[1, 2] == 1 && B[1, 3] == 0
   @test B[2, 1] == 0 && B[2, 2] == 1 && B[2, 3] == 1

   println("PASS")
end =#

function test_sresolution_minres()
   print("sresolution.minres...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x*y + 1, x^2 + 1)
   F = fres(std(I), 3)

   M = minres(F)

   @test length(M) == 2

   M1 = Singular.Matrix(M[1])
   M2 = Singular.Matrix(M[2])

   # check it is a complex
   @test iszero(M1*M2)

   println("PASS")
end

function test_sresolution()
   test_sresolution_constructors()
   test_sresolution_manipulation()
#    test_sresolution_betti()
   test_sresolution_minres()

   println("")
end

