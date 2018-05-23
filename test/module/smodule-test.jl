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

   println("PASS")
end

function test_smodule()
   test_smodule_constructors()
   test_smodule_manipulation()

   println("")
end
