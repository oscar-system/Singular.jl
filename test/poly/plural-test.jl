@testset "plural.constructors" begin

   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

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

   @test symbols(R) == [:x, :y]
   @test Singular.singular_symbols(R) == symbols(R)
end

@testset "plural.printing" begin
   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

   @test string(x) == "x"
   @test string(y) == "y"
   @test string(x^2 + 2*y + 3) == "x^2 + 2*y + 3"
end

@testset "sglpoly.rename" begin
   s = ["x[1]", "x[2]"]
   r, x = polynomial_ring(QQ, s, ordering = :degrevlex)
   R, x = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                      Singular.Matrix(r, [0 x[1]; 0 0]))
   @test String.(symbols(R)) == s
   @test String.(Singular.singular_symbols(R)) == ["x_1", "x_2"]
end

@testset "plural.manipulation" begin
   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

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

   r, (x, y) = polynomial_ring(Fp(5), ["x", "y"], ordering = :degrevlex)
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

@testset "plural.binary_ops" begin

   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))
   @test y*x == x*y + x

   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x+1; 0 0]))

   @test x * QQ(2) == QQ(2) * x
   @test x * 2 == 2 * x
   @test y*x == x*y + x + 1
end

@testset "plural.hash" begin
   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end

@testset "plural.recognition" begin
   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

   S = Singular.create_ring_from_singular_ring(Singular.libSingular.rCopy(R.ptr))
   @test S isa PluralRing{n_Q}
   @test nvars(S) == 2
end

@testset "plural.ideal" begin
   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x; 0 0]))

   S = Singular.create_ring_from_singular_ring(Singular.libSingular.rCopy(R.ptr))

   # zero ideal
   I = Ideal(S)
   @test base_ring(I) isa base_ring_type(I)
   @test parent(I) isa parent_type(I)

   # principal ideal
   I = Ideal(S, gen(S,1))
   @test base_ring(I) isa base_ring_type(I)
   @test parent(I) isa parent_type(I)
end
