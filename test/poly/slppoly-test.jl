@testset "slppoly.constructors" begin

   R, (x, y, z) = FreeAlgebra(QQ, ["x", "y", "z"], 2)

   @test elem_type(R) == slppoly{n_Q}
   @test elem_type(LPPolyRing{n_Q}) == slppoly{n_Q}
   @test parent_type(slppoly{n_Q}) == LPPolyRing{n_Q}
   @test base_ring(R) == QQ

   @test R isa Nemo.AbstractAlgebra.NCRing

   a = R()

   @test base_ring(a) == QQ
   @test parent(a) == R

   @test isa(a, slppoly)

   b = R(123)

   @test isa(b, slppoly)

   c = R(BigInt(123))

   @test isa(c, slppoly)

   d = R(c)

   @test isa(d, slppoly)

   f = R(Nemo.ZZ(123))

   @test isa(f, slppoly)

   g = R(ZZ(123))

   @test isa(g, slppoly)

   @test isgen(x)
   @test isgen(y)
   @test !isgen(R(1))
   @test !isgen(R(0))
   @test !isgen(2x)
   @test !isgen(x + y)
   @test !isgen(x^2)

   @test has_global_ordering(R)
   @test !has_local_ordering(R)
   @test !has_mixed_ordering(R)

   @test symbols(R) == [:x, :y, :z]
end

@testset "slppoly.printing" begin
end

@testset "slppoly.rename" begin
end

@testset "slppoly.manipulation" begin
   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 5)

   @test isone(one(R))
   @test iszero(zero(R))
   @test isunit(R(1))
   @test isunit(R(-1))
   @test isunit(R(2))
   @test !isunit(R(0))
   @test !isunit(x)
   @test isgen(x)
   @test !isgen(R(1))
   @test !isgen(x + 1)
   @test isterm(2x)
   @test !isterm(x + 1)
   @test length(x^2 + 2x + 1) == 3

   @test deepcopy(x + 2) == x + 2

   @test characteristic(R) == 0

   @test nvars(R) == 2
   pol = x^5 - x*y*x + 3*y*x + 2

   @test collect(coefficients(pol)) == [QQ(1), QQ(-1), QQ(3), QQ(2)]
   @test collect(exponent_words(pol)) == [[1,1,1,1,1], [1,2,1], [2,1], Int[]]

   polzip = zip(coefficients(pol), monomials(pol), terms(pol))
   r = R()
   for (c, m, t) in polzip
      r += c*m
      @test t == c*m
   end
   @test pol == r

   B = MPolyBuildCtx(R)
   push_term!(B, QQ(2), [1,2])
   push_term!(B, QQ(3), Int[])
   @test finish(B) == 2*x*y + 3
   @test finish(B) == 0
   push_term!(B, QQ(-1), [1,2,2,1])
   @test finish(B) == -x*y*y*x
   @test finish(B) == 0
   @test_throws Exception push_term!(B, QQ(2), [3,2])
   @test_throws Exception push_term!(B, QQ(2), [1,1,1,1,1,1])

   R, (x, ) = FreeAlgebra(Fp(5), ["x", ], 1)

   @test characteristic(R) == 5

   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 10)
   p = x + y
   q = x

   @test x * QQ(2) == QQ(2) * x
   @test x * 2 == 2 * x

   for i = 1:nvars(R)
      @test gen(R, i) == gens(R)[i]
   end
   @test gen(R, 1) == x

   @test tail(3x^2*y + 2x*y + y + 7) == 2x*y + y + 7
   @test tail(R(1)) == 0
   @test tail(R()) == 0
   @test leading_coefficient(zero(R)) == 0
   @test leading_coefficient(3x^2 + 2x + 1) == 3
   #@test constant_coefficient(x^2*y + 2x + 3) == 3
   #@test constant_coefficient(x^2 + y) == 0
   @test leading_monomial(3x^2 + 2x + 1) == x^2
   @test leading_term(3x^2 + 2x + 1) == 3x^2
   @test trailing_coefficient(3x^2*y + 2x + 7y + 9) == 9
   @test trailing_coefficient(5x) == 5
   @test trailing_coefficient(R(3)) == 3
   @test trailing_coefficient(R()) == 0
end

@testset "slppoly.change_base_ring" begin
end

@testset "slppoly.multivariate_coeff" begin
end

@testset "slppoly.unary_ops" begin
end

@testset "slppoly.binary_ops" begin
   R, (x, y, z, w) = FreeAlgebra(QQ, ["x", "y", "z", "w"], 5)
   @test y*x != x*y
   @test !iszero(y*x + x*y)
   @test !iszero(y*x - x*y)
end

@testset "slppoly.comparison" begin
end

@testset "slppoly.powering" begin
   R, (x, y, z) = FreeAlgebra(Fp(5), ["x", "y", "z"], 10)
   @test length((x + y + z)^5) == 3^5
end

@testset "slppoly.exact_division" begin
end

@testset "slppoly.adhoc_exact_division" begin
end

@testset "slppoly.adhoc_binary_operation" begin
end

@testset "slppoly.divides" begin
end

@testset "slppoly.hash" begin
   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 5)

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end
