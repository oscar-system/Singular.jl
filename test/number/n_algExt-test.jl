@testset "n_algExt.constructors" begin
   F, (a,) = FunctionField(QQ, ["a"])
   K, _ = AlgebraicExtensionField(F, a^2 + 1)
   K1, _ = AlgebraicExtensionField(F, a^2 + 1)
   K2, _ = AlgebraicExtensionField(F, a^2 + 1, cached = false)

   @test K isa Singular.Field
   @test K1 isa Singular.Field
   @test K2 isa Singular.Field
   @test K == K1
   @test K != K2
   @test K1 != K2

   @test_throws ArgumentError AlgebraicExtensionField(F, zero(F))

   F, (a, b) = FunctionField(QQ, ["a", "b"])

   @test_throws ArgumentError AlgebraicExtensionField(F, a^2 + b)
end

@testset "n_algExt.printing" begin
   F, (a,) = FunctionField(QQ, ["a"])
   F, a = AlgebraicExtensionField(F, a^2 + 1)

   @test string(a+1) == "(a+1)"

   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])

   p = (3*a+2)*x+(2*a^2)*y+a*z
   q = (a+1)*x^2
   # singular's printing
   @test string(p) == "(3*a+2)*x-2*y+(a)*z"
   @test string(q) == "(a+1)*x^2"
   # AA's printing
   @test sprint(show, "text/plain", p) == "(3*a + 2)*x - 2*y + a*z"
   @test sprint(show, "text/plain", q) == "(a + 1)*x^2"
end

@testset "n_algExt.manipulation" begin
   F, (a,) = FunctionField(Fp(7), ["a"])
   K, _ = AlgebraicExtensionField(F, a^2 + 1)

   @test isone(one(K))
   @test iszero(zero(K))
   @test isunit(K(1)) && isunit(K(2))
   @test !isunit(K(0))

   @test characteristic(K) == 7
   @test modulus(K) == a^2 + 1

   @test deepcopy(K(2)) == K(2)

   a = gen(K)
   d = Dict{elem_type(K), Int}()
   d[a] = 1
   d[a + 1] = 2
   d[a] = 3
   @test d[a + 1] == 2
   @test d[a] == 3
end

@testset "n_algExt.unary_ops" begin
   F, (a,) = FunctionField(Fp(7), ["a"])
   F, a = AlgebraicExtensionField(F, a^2 + 1)

   @test -F(3) == F(4)
   @test -F() == F()
   @test a^3 == -a
end

@testset "n_algExt.binary_ops" begin
   F, (a,) = FunctionField(QQ, ["a"])
   F, a = AlgebraicExtensionField(F, a^2 + 1)

   x = 2*a + 2
   y = 3*a + 2

   @test x + y == 5*a + 4
   @test x - y == -a
   @test x*y == 10*a - 2
end

@testset "n_algExt.comparison" begin
   F, (a,) = FunctionField(Fp(7), ["a"])
   F, a = AlgebraicExtensionField(F, a^2 + 1)

   x = -2*a - 3

   @test x == deepcopy(x)
   @test isequal(x, x)
end

@testset "n_algExt.powering" begin
   F, (a,) = FunctionField(Fp(7), ["a"])
   F, a = AlgebraicExtensionField(F, a^2 + 1)

   @test (a + 1)^3 == 2*a - 2
   @test -4*(a + 1)^-3 == a + 1

   @test_throws DivideError (a - a)^-2
end

@testset "n_algExt.exact_division" begin
   F, (a,) = FunctionField(QQ, ["a"])
   F, a = AlgebraicExtensionField(F, a^2 + 1)

   x = a^2 - 1
   y = a + 1

   @test inv(x) == 1 // (a^2 - 1)
   @test divexact(x, y) == a - 1
end

@testset "n_algExt.gcd_lcm" begin
   F, (a,) = FunctionField(Fp(7), ["a"])
   F, a = AlgebraicExtensionField(F, a^2 + 1)

   x = a^2 - 1
   y = a + 1

   @test gcd(x, y) == F(1)
   @test gcd(F(0), F(0)) == F(0)
end

@testset "n_algExt.recognition" begin

   F, (a,) = FunctionField(QQ,["a"])
   K, a = AlgebraicExtensionField(F, a^2 + 1)
   R, (x, y, z) = PolynomialRing(K, ["x", "y", "z"])
   S = Singular.create_ring_from_singular_ring(Singular.libSingular.rCopy(R.ptr))

   @test isa(S, Singular.PolyRing{n_algExt})

   (x, y, z) = gens(S)
   a = gen(base_ring(S))

   @test isa(a*x + y, spoly{n_algExt})
   @test iszero((a*x)^2 + x^2)

   @test base_ring(S) == K
   @test base_ring(a*x) == K
   @test parent(a) == K
end
