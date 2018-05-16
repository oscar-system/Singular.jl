function test_svector_constructors()
   print("svector.constructors...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   S1 = parent(v1)

   M = FreeModule(R, 3)
   v2 = M([x + 1, x*y + 1, y])
   S2 = parent(v2)

   @test elem_type(S1) == svector{spoly{n_Q}}
   @test elem_type(FreeMod{spoly{n_Q}}) == svector{spoly{n_Q}}
   @test parent_type(svector{spoly{n_Q}}) == FreeMod{spoly{n_Q}}

   @test base_ring(S1) == R
   @test base_ring(v1) == R

   @test typeof(S1) <: Singular.Module
   @test typeof(S1) <: AbstractAlgebra.Module
   @test typeof(v1) <: AbstractAlgebra.ModuleElem

   @test isa(v1, svector)

   println("PASS")
end

function test_svector_manipulation()
   print("svector.manipulation...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   M = FreeModule(R, 3)

   v = vector(R, x, y, R(2))

   @test rank(M) == 3

   g = gens(M)

   @test length(g) == 3

   @test g[1]*x + g[2]*y + g[3]*2 == v

   @test deepcopy(v) == v

   println("PASS")
end

function test_svector_unary_ops()
   print("svector.unary_ops...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v = vector(R, x, y, R(2))

   @test -(-v) == v

   println("PASS")
end

function test_svector_binary_ops()
   print("svector.binary_ops...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x, y, R(2))
   v2 = vector(R, x^2 + 1, x*y, y + 1)

   @test v1 + v2 == v2 + v1
   @test v1 + v2 - v2 == v1

   println("PASS")
end

function test_svector_adhoc_binary()
   print("svector.adhoc_binary...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x, y, R(2))
   v2 = vector(R, x^2 + 1, x*y, y + 1)

   @test 2*(v1 + v2) == v1*2 + v2*2
   @test QQ(2)*(v1 + v2) == v1*QQ(2) + v2*QQ(2)
   @test x*(v1 + v2) == v1*x + v2*x
   
   println("PASS")
end

function test_svector_comparison()
   print("svector.comparison...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x, y, R(2))
   v2 = vector(R, x^2 + 1, x*y, y + 1)

   @test v1 != v2
   @test v1 == deepcopy(v1)

   println("PASS")
end

function test_svector()
   test_svector_constructors()
   test_svector_manipulation()
   test_svector_unary_ops()
   test_svector_binary_ops()
   test_svector_adhoc_binary()
   test_svector_comparison()

   println("")
end
