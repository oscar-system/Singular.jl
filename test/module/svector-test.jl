@testset "svector.constructors..." begin
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
end

@testset "svector.jet..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   a = vector(R, x^5 + 1, 2x^3 + 3y^2, x^2)

   b = vector(R, R(1), 2*x^3 + 3*y^2, x^2)

   c = jet(a, 3)

   @test b == c
end

@testset "svector.manipulation..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   M = FreeModule(R, 3)

   v = vector(R, x, y, R(2))

   @test rank(M) == 3

   g = gens(M)

   @test length(g) == 3

   @test g[1]*x + g[2]*y + g[3]*2 == v

   @test deepcopy(v) == v
end

@testset "svector.unary_ops..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v = vector(R, x, y, R(2))

   @test -(-v) == v
end

@testset "svector.binary_ops..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x, y, R(2))
   v2 = vector(R, x^2 + 1, x*y, y + 1)

   @test v1 + v2 == v2 + v1
   @test v1 + v2 - v2 == v1
end

@testset "svector.adhoc_binary..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x, y, R(2))
   v2 = vector(R, x^2 + 1, x*y, y + 1)

   @test 2*(v1 + v2) == v1*2 + v2*2
   @test QQ(2)*(v1 + v2) == v1*QQ(2) + v2*QQ(2)
   @test x*(v1 + v2) == v1*x + v2*x
end

@testset "svector.comparison..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x, y, R(2))
   v2 = vector(R, x^2 + 1, x*y, y + 1)

   @test v1 != v2
   @test v1 == deepcopy(v1)
end

@testset "svector.conversion..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   v1 = vector(R, x, y, R(2))

   @test Array(v1) == [x, y, R(2)]
end
