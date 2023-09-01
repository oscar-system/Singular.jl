@testset "sresolution.constructors" begin
   R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])

   I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
   F = fres(std(I), 3)
   S = parent(F)

   @test elem_type(S) == sresolution{spoly{n_Q}}
   @test elem_type(ResolutionSet{spoly{n_Q}}) == sresolution{spoly{n_Q}}
   @test parent_type(sresolution{spoly{n_Q}}) == ResolutionSet{spoly{n_Q}}

   @test base_ring(S) == R
   @test base_ring(F) == R

   @test S isa AbstractAlgebra.Set
   @test F isa AbstractAlgebra.SetElem

   @test isa(F, sresolution)
end

@testset "sresolution.manipulation" begin
   R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])

   I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
   F = fres(std(I), 1)

   @test length(F) == 1

   F = fres(std(I), 0)

   @test length(F) == 3

   F = fres(std(I), 3)

   @test length(F) == 3

   F = fres(std(I), 4)

   @test length(F) == 3

   @test isa(F[1], sideal)
   @test isa(F[2], smodule)
   @test isa(F[3], smodule)

   G = deepcopy(F)

   @test isa(G, sresolution)
   @test parent(F) === parent(G)
end

@testset "sresolution.betti" begin
   R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])

   I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
   F = fres(std(I), 3)
   M = minres(F)

   B = betti(M)

   @test size(B) == (3, 4)

   @test B[1, 1] == 1 && B[1, 2] == 0 && B[1, 3] == 0 && B[1, 4] == 0
   @test B[2, 1] == 0 && B[2, 2] == 5 && B[2, 3] == 5 && B[2, 4] == 0
   @test B[3, 1] == 0 && B[3, 2] == 0 && B[3, 3] == 0 && B[3, 4] == 1
end

@testset "sresolution.minres" begin
   R, (w, x, y, z) = polynomial_ring(QQ, ["w", "x", "y", "z"])

   I = Ideal(R, w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z)
   F = fres(std(I), 3)

   M = minres(F)

   @test length(M) == 3

   M1 = Singular.Matrix(M[1])
   M2 = Singular.Matrix(M[2])

   # check it is a complex
   @test iszero(M1*M2)
end

@testset "sresolution.fres_module" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

   I = Singular.Ideal(R, y*z + z^2, y^2 + x*z, x*y + z^2, z^3, x*z^2, x^2*z)
   M = std(syz(I))

   F = fres(M, 0)

   # We have R^6 <- R^9 <- R^5 <- R^1
   # All references agree that when written as follows:
   # 0 <- M <- R^6 <- R^9 <- R^5 <- R^1 <- 0
   # the length is 3, as numbering starts at index 0 with R^6

   @test length(F) == 3
   @test F[1] isa smodule
   @test F[2] isa smodule

   M1 = Singular.Matrix(F[1])
   M2 = Singular.Matrix(F[2])
   M3 = Singular.Matrix(F[3])

   @test Singular.Matrix(Singular.Module(M1)) == M1

   @test iszero(M1*M2)
   @test iszero(M2*M3)
end

@testset "sresolution.mres_module" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, y )
   v2 = vector(R, x)
   v3 = vector(R, x+y)

   M = Singular.Module(R, v1, v2, v3)
   L,TT = mres_with_map(M,0)
   @test iszero(Singular.Matrix(M)*TT-Singular.Matrix(L[1]))
   @test iszero(Singular.Matrix(L[1])*Singular.Matrix(L[2]))
end
