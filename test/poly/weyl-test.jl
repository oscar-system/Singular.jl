@testset "weyl.constructors" begin

   R, (x, y, z, dx, dy, dz) = WeylAlgebra(QQ, ["x", "y", "z"])

   @test typeof(R) == PluralRing{n_Q}
   @test elem_type(R) == spluralg{n_Q}
   @test elem_type(PluralRing{n_Q}) == spluralg{n_Q}
   @test parent_type(spluralg{n_Q}) == PluralRing{n_Q}
   @test base_ring(R) == QQ

   @test R isa Nemo.AbstractAlgebra.NCRing

   a = R()

   @test base_ring(a) == QQ
   @test parent(a) == R

   @test isa(a, spluralg)

   b = R(123)

   @test isa(b, spluralg)

   c = R(BigInt(123))

   @test isa(c, spluralg)

   d = R(c)

   @test isa(d, spluralg)

   f = R(Nemo.ZZ(123))

   @test isa(f, spluralg)

   g = R(ZZ(123))

   @test isa(g, spluralg)


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

   @test symbols(R) == [:x, :y, :z, :dx, :dy, :dz]
   @test Singular.singular_symbols(R) == symbols(R)

   @test_throws Exception WeylAlgebra(QQ, ["x", "y"], ordering = :neglex)
end

@testset "weyl.printing" begin
   r, (x, y, dx, dy) = WeylAlgebra(QQ, ["x", "y"], ordering = :degrevlex)

   @test string(x) == "x"
   @test string(y) == "y"
   @test string(dx) == "dx"
   @test string(dy) == "dy"
   @test string(x^2*dy + 2*y + 3) == "x^2*dy + 2*y + 3"
end

@testset "weyl.rename" begin
   R, x = WeylAlgebra(QQ, ["x[1]", "x[2]", "x[3]"])
   @test String.(symbols(R)) == ["x[1]", "x[2]", "x[3]", "dx[1]", "dx[2]", "dx[3]"]
   @test String.(Singular.singular_symbols(R)) == ["x_1", "x_2", "x_3", "dx_1", "dx_2", "dx_3"]
end

@testset "weyl.manipulation" begin
   R, (x, dx) = WeylAlgebra(QQ, ["x"])

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
   pol = x^5 + 3x + 2

   @test collect(coefficients(pol)) == [QQ(1), QQ(3), QQ(2)]
   @test collect(exponent_vectors(pol)) == [[5,0], [1,0], [0,0]]

   polzip = zip(coefficients(pol), monomials(pol), terms(pol))
   r = R()
   for (c, m, t) in polzip
      r += c*m
      @test t == c*m
   end
   @test pol == r

   B = MPolyBuildCtx(R)
   push_term!(B, QQ(2), [1,2])
   push_term!(B, QQ(3), [0,0])
   @test finish(B) == 2*x*dx^2 + 3
   @test finish(B) == 0
   push_term!(B, QQ(-1), [10,0])
   @test finish(B) == -x^10
   @test finish(B) == 0
   @test_throws Exception push_term!(B, QQ(2), [0,0,0])

   R, (x, dx) = WeylAlgebra(Fp(5), ["x"; "dx"])

   @test characteristic(R) == 5

   R, (x, y, dx, dy) = WeylAlgebra(QQ, ["x", "y"])
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
   @test leading_monomial(3x^2 + 2x + 1) == x^2
   @test leading_term(3x^2 + 2x + 1) == 3x^2
   @test trailing_coefficient(3x^2*y + 2x + 7y + 9) == 9
   @test trailing_coefficient(5x) == 5
   @test trailing_coefficient(R(3)) == 3
   @test trailing_coefficient(R()) == 0
end

@testset "weyl.binary_ops" begin
   R, (x, y, dx, dy) = WeylAlgebra(QQ, ["x" "y"; "dx" "dy"])
   @test y*x == x*y
   @test dy*dx == dx*dy
   @test dy*y == y*dy + 1
   @test dx*x == x*dx + 1
   @test dx^2*(x^2 + y) == 2 + 4*x*dx + (x^2 + y)*dx^2
end

@testset "weyl.hash" begin
   R, (x, y) = WeylAlgebra(QQ, ["x", "y"])

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end
