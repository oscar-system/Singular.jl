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
end

@testset "sextpoly.printing" begin
end

@testset "sextpoly.rename" begin
end

@testset "sextpoly.manipulation" begin
   R, (x, ) = ExteriorAlgebra(QQ, ["x", ])

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

   @test nvars(R) == 1
   pol = x^5 + 3x + 2

   @test length(collect(coefficients(pol))) == length(pol)
   @test length(collect(exponent_vectors(pol))) == length(pol)

   polzip = zip(coefficients(pol), monomials(pol), terms(pol))
   r = R()
   for (c, m, t) in polzip
      r += c*m
      @test t == c*m
   end
   @test pol == r

   R, (x, ) = ExteriorAlgebra(Fp(5), ["x", ])

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

   #@test isordering_symbolic(R)
   #@test ordering_as_symbol(R) == :degrevlex
   #@test degree(x^2*y^3 + 1, 1) == 2
   #@test degree(x^2*y^3 + 1, y) == 3
   #@test degree(R(), 1) == -1
   #@test degrees(x^2*y^3) == [2, 3]
   #@test vars(x^2 + 3x + 1) == [x]
   #@test var_index(x) == 1 && var_index(y) == 2
   @test tail(3x^2*y + 2x*y + y + 7) == y + 7
   @test tail(R(1)) == 0
   @test tail(R()) == 0
   @test leading_coefficient(zero(R)) == 0
   @test leading_coefficient(3x^2 + 2x + 1) == 2
   #@test constant_coefficient(x^2*y + 2x + 3) == 3
   #@test constant_coefficient(x^2 + y) == 0
   @test leading_monomial(3x^2 + 2x + 1) == x
   @test leading_term(3x^2 + 2x + 1) == 2x
   @test trailing_coefficient(3x^2*y + 2x + 7y + 9) == 9
   @test trailing_coefficient(5x) == 5
   @test trailing_coefficient(R(3)) == 3
   @test trailing_coefficient(R()) == 0
end

@testset "sextpoly.change_base_ring" begin
end

@testset "sextpoly.multivariate_coeff" begin
end

@testset "sextpoly.unary_ops" begin
end

@testset "sextpoly.binary_ops" begin
   R, (x, y, z, w) = ExteriorAlgebra(QQ, ["x", "y", "z", "w"])
   @test y*x == -x*y
   @test y*x*y == 0
   @test 0 == x^2
   @test (x*y + y*z)*(w + x) == x*y*w + y*z*w + y*z*x
end

@testset "sextpoly.comparison" begin
end

@testset "sextpoly.powering" begin
end

@testset "sextpoly.exact_division" begin
end

@testset "sextpoly.adhoc_exact_division" begin
end

@testset "sextpoly.adhoc_binary_operation" begin
end

@testset "sextpoly.divides" begin
end

@testset "sextpoly.hash" begin
   R, (x, y) = ExteriorAlgebra(QQ, ["x", "y"])

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end
