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
   F = fres(std(I), 3)

   println("PASS")
end

function test_sresolution()
   test_sresolution_constructors()
   test_sresolution_manipulation()

   println("")
end

