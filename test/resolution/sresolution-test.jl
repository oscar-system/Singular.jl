@testset "sresolution.constructors..." begin
   R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])

   I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
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
end

@testset "sresolution.manipulation..." begin
   R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])

   I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
   F = fres(std(I), 1)

   @test length(F) == 1

   F = fres(std(I), 0)

   @test length(F) == 3

   F = fres(std(I), 3)

   @test length(F) == 3

   F = fres(std(I), 4)

   @test length(F) == 3

   @test isa(F[1], smodule)
   @test isa(F[2], smodule)
   @test isa(F[3], smodule)

   G = deepcopy(F)

   @test isa(G, sresolution)
   @test parent(F) === parent(G)
end

@testset "sresolution.betti..." begin
   R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])

   I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
   F = fres(std(I), 3)
   M = minres(F)

   B = betti(M)

   @test size(B) == (3, 4)

   @test B[1, 1] == 1 && B[1, 2] == 0 && B[1, 3] == 0 && B[1, 4] == 0
   @test B[2, 1] == 0 && B[2, 2] == 5 && B[2, 3] == 5 && B[2, 4] == 0
   @test B[3, 1] == 0 && B[3, 2] == 0 && B[3, 3] == 0 && B[3, 4] == 1
end

@testset "sresolution.minres..." begin
   R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])

   I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
   F = fres(std(I), 3)

   M = minres(F)

   @test length(M) == 3

   M1 = Singular.Matrix(M[1])
   M2 = Singular.Matrix(M[2])

   # check it is a complex
   @test iszero(M1*M2)
end

