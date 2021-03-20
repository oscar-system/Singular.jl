@testset "n_transExt.constructors" begin
   F, _ = FunctionField(QQ, ["a", "b", "c"])
   F1, _ = FunctionField(QQ, ["a", "b", "c"])
   F2, _ = FunctionField(QQ, ["a", "b", "c"], cached = false)

   @test F isa Singular.Field
   @test F1 isa Singular.Field
   @test F2 isa Singular.Field
   @test F == F1
   @test F != F2
   @test F1 != F2

   @test_throws ArgumentError FunctionField(QQ, String[])
   @test_throws ArgumentError FunctionField(QQ, ["", "b"])
   @test_throws ArgumentError FunctionField(QQ, ["a", ""])
   @test_throws ArgumentError FunctionField(QQ, [""])
   @test_throws ArgumentError FunctionField(QQ, ["a", "b", "a"])
end

@testset "n_transExt.printing" begin
   F, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])

   t = 3*a*b + 2*c

   @test sprint(show, "text/plain", t) == "3*a*b + 2*c"
   @test string(t) == "3*a*b + 2*c"

   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])

   p = (3*a*b+2*c)*x+(c^2)*y+b*z
   q = (a+b//c)*x^2

   @test sprint(show, "text/plain", p) == "(3*a*b + 2*c)*x + c^2*y + b*z"
   @test string(p) == "(3*a*b + 2*c)*x + c^2*y + b*z"

   @test sprint(show, "text/plain", q) == "(a*c + b)//c*x^2"
   @test string(q) == "(a*c + b)//c*x^2"

   @test sprint(show, "text/plain", R(3)) == "3"
   @test string(R(3)) == "3"

   @test sprint(show, "text/plain", R(big(3)^80)) == string(big(3)^80)
   @test string(R(big(3)^80)) == string(big(3)^80)
end

@testset "n_transExt.manipulation" begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])
   x = (a*b+c^2) // (a+b)

   R, (x1, x2, x3) = PolynomialRing(Fp(5), ["x1", "x2", "x3"])
   p = x1*x2 + x3^2

   @test numerator(x) == a*b+c^2
   @test denominator(x) == a+b

   @test isone(one(F))
   @test iszero(zero(F))
   @test isunit(F(1)) && isunit(F(2))
   @test !isunit(F(0))

   @test characteristic(F) == 5
   @test transcendence_degree(F) == 3

   @test deepcopy(F(2)) == F(2)
   @test n_transExt_to_spoly(x, parent_ring = R) == p
end

@testset "n_transExt.unary_ops" begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])
   x = 3*a + b

   @test -F(3) == F(2)
   @test -F() == F()
   @test -x == 2*a - b
end

@testset "n_transExt.binary_ops" begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = 2*a + 2*c
   y = 3*b + 2*c

   @test x + y == 2*a + 3*b - c
   @test x - y == 2*a + 2*b
   @test x*y == a*b - a*c + b*c - c^2
end

@testset "n_transExt.comparison" begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2*b*c + b^2*c + 5*a

   @test x == deepcopy(x)
   @test isequal(x, x)
end

@testset "n_transExt.powering" begin
   F, (a, b, c) = FunctionField(Fp(3), ["a", "b", "c"])

   x = a*b*c + 1
   @test x^3 == a^3*b^3*c^3 + 1

   @test_throws DomainError x^-rand(1:99)
end

@testset "n_transExt.exact_division" begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2 - 1
   y = a + 1

   @test inv(x) == 1 // (a^2 - 1)
   @test divexact(x, y) == a-1
end

@testset "n_transExt.gcd_lcm" begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2 - 1
   y = a + 1

   @test gcd(x, y) == F(1)
   @test gcd(F(0), F(0)) == F(0)
end

@testset "n_transExt.Polynomials" begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])
   R, (x, ) = PolynomialRing(F, ["x"])

   f = (1 + 2x)*a^2 + 3x*b*c + (x + 1)

   g = f^5

   @test g == (2*a^10-2*b^5*c^5+1)*x^5+(a^10+1)

end
