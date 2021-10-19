@testset "sextpoly.constructors" begin

   R, (x, y, z) = ExteriorAlgebra(QQ, ["x", "y", "z"])

   @test elem_type(R) == sextpoly{n_Q}
   @test elem_type(ExtPolyRing{n_Q}) == sextpoly{n_Q}
   @test parent_type(sextpoly{n_Q}) == ExtPolyRing{n_Q}
   @test base_ring(R) == QQ

   @test R isa Nemo.AbstractAlgebra.NCRing

   a = R()

   @test base_ring(a) == QQ
   @test parent(a) == R

   @test isa(a, sextpoly)

   b = R(123)

   @test isa(b, sextpoly)

   c = R(BigInt(123))

   @test isa(c, sextpoly)

   d = R(c)

   @test isa(d, sextpoly)

   f = R(Nemo.ZZ(123))

   @test isa(f, sextpoly)

   g = R(ZZ(123))

   @test isa(g, sextpoly)


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

   @test length(symbols(R)) == 3
   @test symbols(R) == [:x, :y, :z]
   @test Singular.singular_symbols(R) == symbols(R)
end

@testset "sweylpoly.printing" begin
   r, (dx, dy, dz) = ExteriorAlgebra(QQ, ["dx", "dy", "dz"], ordering = :degrevlex)

   @test string(dx) == "dx"
   @test string(dy) == "dy"
   @test string(dz) == "dz"
   @test string(2*dx*dy*dz - 1) == "2*dx*dy*dz - 1"
end

@testset "sweylpoly.rename" begin
   s = ["x[1]", "x[2]", "x[3]"]
   R, x = ExteriorAlgebra(QQ, s)
   @test String.(symbols(R)) == s
   @test String.(Singular.singular_symbols(R)) == ["x_1", "x_2", "x_3"]
end

@testset "sextpoly.manipulation" begin
   R, (x, y) = ExteriorAlgebra(QQ, ["x", "y"])

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
   @test length(x^2 + 2x + 1) == 2

   @test deepcopy(x + 2) == x + 2

   @test characteristic(R) == 0

   @test nvars(R) == 2
   pol = x^5 + 3x + 2

   @test collect(coefficients(pol)) == [QQ(3), QQ(2)]
   @test collect(exponent_vectors(pol)) == [[1, 0], [0,0]]

   polzip = zip(coefficients(pol), monomials(pol), terms(pol))
   r = R()
   for (c, m, t) in polzip
      r += c*m
      @test t == c*m
   end
   @test pol == r

   B = MPolyBuildCtx(R)
   push_term!(B, QQ(3), [0,0])
   push_term!(B, QQ(2), [1,0])
   @test finish(B) == 2*x + 3
   @test finish(B) == 0
   push_term!(B, QQ(-1), [10,0])
   @test finish(B) == 0
   @test finish(B) == 0
   push_term!(B, QQ(1), [1,0])
   push_term!(B, QQ(-1), [0,1])
   @test finish(B) == x - y
   @test_throws Exception push_term!(B, QQ(2), [0,0,0])

   R, (x, ) = ExteriorAlgebra(Fp(5), ["x", "y", "z", "w"])

   @test characteristic(R) == 5

   R, (x, y) = ExteriorAlgebra(QQ, ["x", "y"])
   p = x + y
   q = x

   @test x * QQ(2) == QQ(2) * x
   @test x * 2 == 2 * x

   for i = 1:nvars(R)
      @test gen(R, i) == gens(R)[i]
   end
   @test gen(R, 1) == x

   @test tail(3x^2*y + 2x*y + y + 7) == y + 7
   @test tail(R(1)) == 0
   @test tail(R()) == 0
   @test leading_coefficient(zero(R)) == 0
   @test leading_coefficient(3x^2 + 2x + 1) == 2
   @test constant_coefficient(x^2*y + 2x + 3) == 3
   @test constant_coefficient(x^2 + y) == 0
   @test constant_coefficient(zero(R)) == 0
   @test leading_monomial(3x^2 + 2x + 1) == x
   @test leading_term(3x^2 + 2x + 1) == 2x
   @test trailing_coefficient(3x^2*y + 2x + 7y + 9) == 9
   @test trailing_coefficient(5x) == 5
   @test trailing_coefficient(R(3)) == 3
   @test trailing_coefficient(R()) == 0
end

@testset "sextpoly.binary_ops" begin
   R, (x, y, z, w) = ExteriorAlgebra(QQ, ["x", "y", "z", "w"])
   @test y*x == -x*y
   @test y*x*y == 0
   @test 0 == x^2
   @test (x*y + y*z)*(w + x) == x*y*w + y*z*w + y*z*x
end

@testset "sextpoly.hash" begin
   R, (x, y) = ExteriorAlgebra(QQ, ["x", "y"])

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end
