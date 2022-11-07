@testset "lp.constructors" begin

   R, (x, y, z) = FreeAlgebra(QQ, ["x", "y", "z"], 2)
   R1, _ = FreeAlgebra(QQ, [:x, :y, :z], 2)
   R2, _ = FreeAlgebra(QQ, [:x, :y, :z], 3)
   @test R == R1
   @test R != R2

   @test elem_type(R) == slpalg{n_Q}
   @test elem_type(LPRing{n_Q}) == slpalg{n_Q}
   @test parent_type(slpalg{n_Q}) == LPRing{n_Q}
   @test base_ring(R) == QQ

   @test R isa Nemo.AbstractAlgebra.NCRing

   a = R()

   @test base_ring(a) == QQ
   @test parent(a) == R

   @test isa(a, slpalg)

   b = R(123)

   @test isa(b, slpalg)

   c = R(BigInt(123))

   @test isa(c, slpalg)

   d = R(c)

   @test isa(d, slpalg)

   f = R(Nemo.ZZ(123))

   @test isa(f, slpalg)

   g = R(ZZ(123))

   @test isa(g, slpalg)

   @test is_gen(x)
   @test is_gen(y)
   @test !is_gen(R(1))
   @test !is_gen(R(0))
   @test !is_gen(2x)
   @test !is_gen(x + y)
   @test !is_gen(x^2)

   @test has_global_ordering(R)
   @test !has_local_ordering(R)
   @test !has_mixed_ordering(R)

   @test symbols(R) == [:x, :y, :z]
   @test Singular.singular_symbols(R) == symbols(R)
end

@testset "lp.printing" begin
   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 4, ordering = :degrevlex)

   @test string(x) == "x"
   @test string(y) == "y"
   @test string(zero(R)) == "0"
   @test string(one(R)) == "1"
   @test string(x^2 + 2*y + 3) == "x^2 + 2*y + 3"
   @test string(2*y^2*x^2 + y*x^2*y + x + 3) == "2*y^2*x^2 + y*x^2*y + x + 3"
end

@testset "lp.rename" begin
   for n in 1:5, d in 1:5
      s = ["x[$i]" for i in 1:n]
      R, x = FreeAlgebra(QQ, s, d)
      @test String.(symbols(R)) == s
      @test String.(Singular.singular_symbols(R)) == ["x_$i" for i in 1:n]
   end
end

@testset "lp.manipulation" begin
   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 5)

   @test isone(one(R))
   @test iszero(zero(R))
   @test is_unit(R(1))
   @test is_unit(R(-1))
   @test is_unit(R(2))
   @test !is_unit(R(0))
   @test !is_unit(x)
   @test is_gen(x)
   @test !is_gen(R(1))
   @test !is_gen(x + 1)
   @test is_constant(R(0))
   @test is_constant(R(1))
   @test !is_constant(x)
   @test !is_constant(x + 1)
   @test is_term(2x)
   @test !is_term(x + 1)
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
   @test constant_coefficient(x^2*y + 2x + 3) == 3
   @test constant_coefficient(x^2 + y) == 0
   @test constant_coefficient(zero(R)) == 0
   @test leading_monomial(3x^2 + 2x + 1) == x^2
   @test leading_term(3x^2 + 2x + 1) == 3x^2
   @test trailing_coefficient(3x^2*y + 2x + 7y + 9) == 9
   @test trailing_coefficient(5x) == 5
   @test trailing_coefficient(R(3)) == 3
   @test trailing_coefficient(R()) == 0
end

@testset "lp.binary_ops" begin
   R, (x, y, z, w) = FreeAlgebra(QQ, ["x", "y", "z", "w"], 5)
   @test y*x != x*y
   @test !iszero(y*x + x*y)
   @test !iszero(y*x - x*y)
end

@testset "lp.powering" begin
   R, (x, y, z) = FreeAlgebra(Fp(5), ["x", "y", "z"], 10)
   @test length((x + y + z)^5) == 3^5
end

@testset "lp.hash" begin
   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 5)

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end

@testset "lp.recognition" begin
   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 5)
   S = Singular.create_ring_from_singular_ring(Singular.libSingular.rCopy(R.ptr))
   @test S isa LPRing{n_Q}
   @test degree_bound(S) == 5
   @test nvars(S) == 2
end

