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

   println("PASS")
end

function test_svector()
   test_svector_constructors()
   test_svector_manipulation()

   println("")
end

