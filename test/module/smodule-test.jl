function test_smodule_constructors()
   print("smodule.constructors...")

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

   println("PASS")
end

function test_smodule_jet()
   print("smodule.jet...")

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

   println("PASS")
end

function test_smodule_local()
   print("smodule.local...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:negdegrevlex)

   v1 = vector(R, x, y^2)
   v2 = vector(R, y - x, y - y^2)
   v3 = v1 + v2

   w1 = vector(R, y, y)
   w2 = vector(R, x - y, y^2 - y)

   M = Singular.Module(R, v1, v2, v3)
   MM = Singular.minimal_generating_set(M)

   @test MM[1] == w1 && MM[2] == w2
   println("PASS")
end

function test_smodule_manipulation()
   print("smodule.manipulation...")

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

   println("PASS")
end

#= function test_smodule_slimgb()
   print("smodule.slimgb...")

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

   println("PASS")
end =#

function test_smodule_std()
   print("smodule.std...")

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

   println("PASS")
end

function test_smodule_syz()
   print("smodule.syz...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, (x + 1)*y, (x*y + 1)*y, y)
   v2 = vector(R, (x + 1)*x, (x*y + 1)*x, x)

   M = Singular.Module(R, v1, v2)

   Z = syz(M)

   @test ngens(Z) == 1

   @test Z[1] == vector(R, x, -y)

   println("PASS")
end

function test_smodule_sres()
   print("smodule.sres...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)

   M = std(Singular.Module(R, v1, v2))

   F = sres(M, 0)

   @test length(F) == 1

   M1 = Singular.Matrix(M)
   M2 = Singular.Matrix(F[1])

   @test iszero(M1*M2)

   F = sres(M, 1)

   @test length(F) == 1

   M1 = Singular.Matrix(M)
   M2 = Singular.Matrix(F[1])

   @test iszero(M1*M2)

   F = sres(M, 3)

   @test length(F) == 1

   M1 = Singular.Matrix(M)
   M2 = Singular.Matrix(F[1])

   @test iszero(M1*M2)

   println("PASS")
end

function test_smodule()
   test_smodule_constructors()
   test_smodule_jet()
   test_smodule_local()
   test_smodule_manipulation()
#    test_smodule_slimgb()
   test_smodule_std()
   test_smodule_syz()
   test_smodule_sres()

   println("")
end
