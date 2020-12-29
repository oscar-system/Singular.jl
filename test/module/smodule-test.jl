@testset "smodule.constructors..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)

   M = Singular.Module(R, v1, v2)

   S = parent(M)

   @test elem_type(S) == smodule{spoly{n_Q}}
   @test elem_type(ModuleClass{spoly{n_Q}}) == smodule{spoly{n_Q}}
   @test parent_type(smodule{spoly{n_Q}}) == ModuleClass{spoly{n_Q}}

   @test base_ring(S) == R
   @test base_ring(M) == R

   @test typeof(S) <: AbstractAlgebra.Set
   @test typeof(M) <: AbstractAlgebra.Module

   @test isa(M, smodule)
end

@testset "smodule.jet..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^5 + 1, 2x^3 + 3y^2, x^2)

   w1 = vector(R, x+1, x*y+1, y)
   w2 = vector(R, R(1), 2*x^3 + 3*y^2, x^2)

   M = Singular.Module(R, v1, v2)

   N = jet(M,3)

   P = Singular.Module(R, w1, w2)

   @test P[1] == N[1]
   @test P[2] == N[2]
end

@testset "smodule.local..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:negdegrevlex)

   v1 = vector(R, x, y^2)
   v2 = vector(R, y - x, y - y^2)
   v3 = v1 + v2

   w1 = vector(R, y, y)
   w2 = vector(R, x - y, y^2 - y)

   M = Singular.Module(R, v1, v2, v3)
   MM = Singular.minimal_generating_set(M)

   @test MM[1] == w1 && MM[2] == w2
end

@testset "smodule.manipulation..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)

   M = Singular.Module(R, v1, v2)

   @test rank(M) == 3

   @test ngens(M) == 2

   @test M[1] == v1
   @test M[2] == v2

   N = deepcopy(M)

   @test ngens(N) == 2

   @test N[1] == v1
   @test N[2] == v2
end

#= @testset "smodule.slimgb..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)
   v3 = x*v1 + y*v2 + vector(R, x, y + 1, y^2)

   M = Singular.Module(R, v1, v2, v3)

   G = slimgb(M; complete_reduction=true)

   @test ngens(G) == 3

   @test G[1] == vector(R, x, y + 1, y^2)
   @test G[2] == vector(R, x + 1, x*y + 1, y)
   @test G[3] == vector(R, x^2 + 1, 2*x + 3*y, x)

   @test G.isGB == true

   # Simply test the interfacing works in this case
   G2 = slimgb(M)

   @test G2.isGB == true
end =#

@testset "smodule.std..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)
   v3 = x*v1 + y*v2 + vector(R, x, y + 1, y^2)

   M = Singular.Module(R, v1, v2, v3)

   G = std(M; complete_reduction=true)

   @test ngens(G) == 3

   @test G[1] == vector(R, x, y + 1, y^2)
   @test G[2] == vector(R, x + 1, x*y + 1, y)
   @test G[3] == vector(R, x^2 + 1, 2*x + 3*y, x)

   @test G.isGB == true

   # Simply test the interfacing works in this case
   G2 = std(M)

   @test G2.isGB == true
end

@testset "smodule.syz..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, (x + 1)*y, (x*y + 1)*y, y)
   v2 = vector(R, (x + 1)*x, (x*y + 1)*x, x)

   M = Singular.Module(R, v1, v2)

   Z = syz(M)

   @test ngens(Z) == 1

   @test Z[1] == vector(R, x, -y)
end

@testset "smodule.modulo..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x)
   v2 = vector(R, y)

   A = Singular.Module(R, v1, v2)
   B = Singular.Module(R, v1)

   M = modulo(A,B)

   @test ngens(M) == 2

   @test M[1] == vector(R, R(1), R(0))
   @test M[2] == vector(R, R(0), x)
end

@testset "smodule.lift..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x)
   v2 = vector(R, y)

   A = Singular.Module(R, v1, v2)
   B = Singular.Module(R, v1)

   M,r = lift(A,B)

   @test M[1] == vector(R,R(1),R(0))
   @test iszero(r[1])
end

@testset "smodule.eliminate..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x)
   v2 = vector(R, y)

   A = Singular.Module(R, v1, v2)

   M = eliminate(A,x)

   @test M[1] == v2
end

@testset "smodule.sres..." begin
   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

   I = Singular.Ideal(R, y*z + z^2, y^2 + x*z, x*y + z^2, z^3, x*z^2, x^2*z)
   M = std(syz(I))

   F = sres(M, 0)

   # We have R^6 <- R^9 <- R^5 <- R^1
   # All references agree that when written as follows:
   # 0 <- M <- R^6 <- R^9 <- R^5 <- R^1 <- 0
   # the length is 3, as numbering starts at index 0 with R^6

   @test length(F) == 3

   M1 = Singular.Matrix(F[1])
   M2 = Singular.Matrix(F[2])
   M3 = Singular.Matrix(F[3])

   @test Singular.Matrix(Singular.Module(M1)) == M1

   @test iszero(M1*M2)
   @test iszero(M2*M3)

   F = sres(M, 1)

   @test length(F) >= 1 # Singular can return more than asked for

   M1 = Singular.Matrix(F[1])
   M2 = Singular.Matrix(F[2])

   @test iszero(M1*M2)

   F = sres(M, 2)

   @test length(F) >= 2 # Singular can return more than asked for

   M1 = Singular.Matrix(F[1])
   M2 = Singular.Matrix(F[2])
   M3 = Singular.Matrix(F[3])

   @test iszero(M1*M2)
   @test iszero(M2*M3)

   F = sres(M, 4)

   @test length(F) == 3
end

