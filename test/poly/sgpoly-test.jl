@testset "sgpoly.constructors" begin

   r, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

   @test elem_type(R) == sgpoly{n_Q}
   @test elem_type(GPolyRing{n_Q}) == sgpoly{n_Q}
   @test parent_type(sgpoly{n_Q}) == GPolyRing{n_Q}
   @test base_ring(R) == QQ

   @test R isa Nemo.AbstractAlgebra.NCRing

   a = R()

   @test base_ring(a) == QQ
   @test parent(a) == R

   @test isa(a, sgpoly)

   b = R(123)

   @test isa(b, sgpoly)

   c = R(BigInt(123))

   @test isa(c, sgpoly)

   d = R(c)

   @test isa(d, sgpoly)

   f = R(Nemo.ZZ(123))

   @test isa(f, sgpoly)

   g = R(ZZ(123))

   @test isa(g, sgpoly)


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

   @test symbols(R) == [:x, :y]
   @test Singular.singular_symbols(R) == symbols(R)
end

@testset "sgpoly.printing" begin
   r, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

   @test string(x) == "x"
   @test string(y) == "y"
   @test string(x^2 + 2*y + 3) == "x^2 + 2*y + 3"
end

@testset "sglpoly.rename" begin
   s = ["x[1]", "x[2]"]
   r, x = PolynomialRing(QQ, s, ordering = :degrevlex)
   R, x = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                      Singular.Matrix(r, [0 x[1]; 0 0]))
   @test String.(symbols(R)) == s
   @test String.(Singular.singular_symbols(R)) == ["x_1", "x_2"]
end

@testset "sgpoly.manipulation" begin
   r, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

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

   r, (x, y) = PolynomialRing(Fp(5), ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

   @test characteristic(R) == 5

   F = base_ring(R)

   B = MPolyBuildCtx(R)
   push_term!(B, F(2), [1,2])
   push_term!(B, F(3), [0,0])
   @test finish(B) == 2*x*y^2 + 3
   @test finish(B) == 0
   push_term!(B, F(-1), [10,0])
   @test finish(B) == -x^10
   @test finish(B) == 0
   @test_throws Exception push_term!(B, F(2), [0,0,0])

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
   @test trailing_coefficient(4x) == 4
   @test trailing_coefficient(R(3)) == 3
   @test trailing_coefficient(R()) == 0
end

@testset "sgpoly.binary_ops" begin

   r, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))
   @test y*x == x*y + x

   r, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x+1; 0 0]))

   @test x * QQ(2) == QQ(2) * x
   @test x * 2 == 2 * x
   @test y*x == x*y + x + 1
end

@testset "sgpoly.hash" begin
   r, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end
