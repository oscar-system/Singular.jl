@testset "n_GF.constructors..." begin
   R, x = FiniteField(7, 2, "x")

   @test elem_type(R) == n_GF
   @test elem_type(N_GField) == n_GF
   @test parent_type(n_GF) == N_GField
   @test base_ring(R) == Union{}
   @test isa(x, n_GF)

   @test R isa Nemo.Field

   a = R()

   @test base_ring(a) == Union{}
   @test parent(a) == R

   @test isa(a, n_GF)

   b = R(3)

   @test isa(b, n_GF)

   c = R(BigInt(3))

   @test isa(c, n_GF)

   # segfault
   # d = R(ZZ(3))

   # @test isa(d, n_GF);

   f = R(c)

   @test isa(f, n_GF)

   f = R(Nemo.ZZ(123))

   @test isa(f, n_GF)

   @test_throws DomainError FiniteField(257, 1, "a")
   @test_throws DomainError FiniteField(2, 16, "a")
   @test_throws DomainError FiniteField(2, 64, "a") # errors even if 2^64 == 0
end

@testset "n_GF.printing..." begin
   R, x = FiniteField(5, 2, "x")

   @test string(zero(R)) == "0"
   @test string(one(R)) == "1"
   @test string(x) == "x"
   @test string(x^2) == "x^2"
   @test string(x^11) == "x^11"

   @test sprint(show, "text/plain", zero(R)) == "0"
   @test sprint(show, "text/plain", one(R)) == "1"
   @test sprint(show, "text/plain", x) == "x"
   @test sprint(show, "text/plain", x^2) == "x^2"
   @test sprint(show, "text/plain", x^11) == "x^11"
end

@testset "n_GF.manipulation..." begin
   R, x = FiniteField(5, 2, "x")

   @test isone(one(R))
   @test iszero(zero(R))
   @test isunit(R(1)) && isunit(R(2))
   @test !isunit(R(0))

   @test gen(R) == x

   @test characteristic(R) == 5
   @test degree(R) == 2

   @test deepcopy(R(2)) == R(2)
end

@testset "n_GF.unary_ops..." begin
   R, x = FiniteField(5, 2, "x")

   a = 3x + 1

   @test -R(3) == R(2)
   @test -R() == R()
   @test -a == 2x - 1
end

@testset "n_GF.binary_ops..." begin
   R, x = FiniteField(5, 2, "x")

   a = 2x + 1
   b = 3x + 2

   @test a + b == R(3)
   @test a - b == -x - 1
   @test a*b == x^19
end

@testset "n_GF.adhoc_binary..." begin
   R, x = FiniteField(5, 2, "x")

   @test R(2) + 3 == R(0)
   @test 2 + R(3) == R(0)
   @test R(2) - 3 == R(4)
   @test 2 - R(3) == R(4)
   @test 2*R(3) == R(1)
   @test R(2)*3 == R(1)
end

@testset "n_GF.comparison..." begin
   R, x = FiniteField(5, 2, "x")

   a = 2x + 3

   @test a == deepcopy(a)
   @test isequal(a, a)
end

@testset "n_GF.ad_hoc_comparison..." begin
   R, x = FiniteField(5, 2, "x")

   @test R(2) == 2
   @test x != 2
   @test 2 == R(2)
   @test 2 != x
   @test isequal(R(2), 2)
   @test isequal(2, R(2))
end

@testset "n_GF.powering..." begin
   R, x = FiniteField(5, 2, "x")

   @test (x + 1)^2 == x^20
   @test_throws DomainError x^-rand(1:99)
end

@testset "n_GF.exact_division..." begin
   R, x = FiniteField(5, 2, "x")

   a = x + 1
   b = 2x + 1

   @test inv(a) == x^2
   @test divexact(a, b) == x^14
end

@testset "n_GF.gcd_lcm..." begin
   R, x = FiniteField(5, 2, "x")

   a = 2x + 1
   b = 3x + 2

   @test gcd(a, b) == R(1)
   @test gcd(R(0), R(0)) == R(0)
end

@testset "n_GF.Polynomials..." begin
   R, x = FiniteField(5, 2, "x")
   S, y = Nemo.PolynomialRing(R, "y")

   f = (1 + 2x)*y^2 + 3x*y + (x + 1)

   g = f^2

   @test g == x^16*y^4+x^9*y^3+x^5*y^2+x^23*y+x^20
end

@testset "n_GF.rand..." begin
   R, x = FiniteField(17, 3, "x")

   @inferred rand(R)
   @test rand(R) isa n_GF
   @test rand(R, 2) isa Vector{n_GF}
   @test rand(R, 2, 3) isa Matrix{n_GF}

   Random.seed!(rng, 0)
   t = rand(rng, R)
   @test t isa n_GF
   s2 = rand(rng, R, 2)
   @test s2 isa Vector{n_GF}
   @test length(s2) == 2
   s3 = rand(rng, R, 2, 3)
   @test s3 isa Matrix{n_GF}
   @test size(s3) == (2, 3)

   Random.seed!(rng, 0)
   @test t == rand(rng, R)
   @test s2 == rand(rng, R, 2)
   @test s3 == rand(rng, R, 2, 3)
end
