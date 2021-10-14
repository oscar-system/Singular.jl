@testset "sweylpoly.constructors" begin

   R, (x, y, z, dx, dy, dz) = WeylAlgebra(QQ, ["x", "y", "z"])

   @test elem_type(R) == sweylpoly{n_Q}
   @test elem_type(WeylPolyRing{n_Q}) == sweylpoly{n_Q}
   @test parent_type(sweylpoly{n_Q}) == WeylPolyRing{n_Q}
   @test base_ring(R) == QQ

   @test R isa Nemo.AbstractAlgebra.NCRing

   a = R()

   @test base_ring(a) == QQ
   @test parent(a) == R

   @test isa(a, sweylpoly)

   b = R(123)

   @test isa(b, sweylpoly)

   c = R(BigInt(123))

   @test isa(c, sweylpoly)

   d = R(c)

   @test isa(d, sweylpoly)

   f = R(Nemo.ZZ(123))

   @test isa(f, sweylpoly)

   g = R(ZZ(123))

   @test isa(g, sweylpoly)


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
end

@testset "sweylpoly.printing" begin
end

@testset "sweylpoly.rename" begin
end

@testset "sweylpoly.manipulation" begin
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

   @test length(collect(coefficients(pol))) == length(pol)
   @test length(collect(exponent_vectors(pol))) == length(pol)

   polzip = zip(coefficients(pol), monomials(pol), terms(pol))
   r = R()
   for (c, m, t) in polzip
      r += c*m
      @test t == c*m
   end
   @test pol == r

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

   #@test isordering_symbolic(R)
   #@test ordering_as_symbol(R) == :degrevlex
   #@test degree(x^2*y^3 + 1, 1) == 2
   #@test degree(x^2*y^3 + 1, y) == 3
   #@test degree(R(), 1) == -1
   #@test degrees(x^2*y^3) == [2, 3]
   #@test vars(x^2 + 3x + 1) == [x]
   #@test var_index(x) == 1 && var_index(y) == 2
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

@testset "sweylpoly.change_base_ring" begin
end

@testset "sweylpoly.multivariate_coeff" begin
end

@testset "sweylpoly.unary_ops" begin
end

@testset "sweylpoly.binary_ops" begin
   R, (x, y, dx, dy) = WeylAlgebra(QQ, ["x" "y"; "dx" "dy"])
   @test y*x == x*y
   @test dy*dx == dx*dy
   @test dy*y == y*dy + 1
   @test dx*x == x*dx + 1
   @test dx^2*(x^2 + y) == 2 + 4*x*dx + (x^2 + y)*dx^2
end

@testset "sweylpoly.comparison" begin
end

@testset "sweylpoly.powering" begin
end

@testset "sweylpoly.exact_division" begin
end

@testset "sweylpoly.adhoc_exact_division" begin
end

@testset "sweylpoly.adhoc_binary_operation" begin
end

@testset "sweylpoly.divides" begin
end

@testset "sweylpoly.hash" begin
   R, (x, y) = WeylAlgebra(QQ, ["x", "y"])

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end
