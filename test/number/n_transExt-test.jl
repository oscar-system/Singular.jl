@testset "n_transExt.constructors..." begin
   F, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])

   @test elem_type(F) == n_transExt
   @test elem_type(N_FField) == n_transExt
   @test parent_type(n_transExt) == N_FField
   @test base_ring(F) == QQ

   @test F isa Singular.Field

   a = F()

   @test base_ring(a) == Union{}
   @test parent(a) == F

   @test isa(a, n_transExt)

   b = F(3)

   @test isa(b, n_transExt)

   c = F(BigInt(3))

   @test isa(c, n_transExt)
   @test b == c

   f = F(c)

   @test isa(f, n_transExt)
   @test b == f

   f = F(ZZ(3))

   @test isa(f, n_transExt)
   @test b == f

   f = F(Nemo.ZZ(3))

   @test isa(f, n_transExt)
   @test b == f

   @test_throws ArgumentError FunctionField(QQ, String[])
   @test_throws ArgumentError FunctionField(QQ, ["", "b"])
   @test_throws ArgumentError FunctionField(QQ, ["a", ""])
   @test_throws ArgumentError FunctionField(QQ, [""])
   @test_throws ArgumentError FunctionField(QQ, ["a", "b", "a"])
end

@testset "n_transExt.printing..." begin
   F, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])

   @test string(3*a*b + 2*c) == "(3*a*b+2*c)"

   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])

   p = (3*a*b+2*c)*x+(c^2)*y+b*z
   q = (a+b//c)*x^2
   # singular's printing
   @test string(p) == "(3*a*b+2*c)*x+(c^2)*y+(b)*z"
   @test string(q) == "(a*c+b)/(c)*x^2"
   # AA's printing
   @test sprint(show, "text/plain", p) == "(3*a*b + 2*c)*x + c^2*y + b*z"
   @test sprint(show, "text/plain", q) == "(a*c + b)//c*x^2"
end

@testset "n_transExt.manipulation..." begin
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

@testset "n_transExt.unary_ops..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])
   x = 3*a + b

   @test -F(3) == F(2)
   @test -F() == F()
   @test -x == 2*a - b
end

@testset "n_transExt.binary_ops..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = 2*a + 2*c
   y = 3*b + 2*c

   @test x + y == 2*a + 3*b - c
   @test x - y == 2*a + 2*b
   @test x*y == a*b - a*c + b*c - c^2
end

@testset "n_transExt.comparison..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2*b*c + b^2*c + 5*a

   @test x == deepcopy(x)
   @test isequal(x, x)
end

@testset "n_transExt.powering..." begin
   F, (a, b, c) = FunctionField(Fp(3), ["a", "b", "c"])

   x = a*b*c + 1
   @test x^3 == a^3*b^3*c^3 + 1

   @test_throws DomainError x^-rand(1:99)
end

@testset "n_transExt.exact_division..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2 - 1
   y = a + 1

   @test inv(x) == 1 // (a^2 - 1)
   @test divexact(x, y) == a-1
end

@testset "n_transExt.gcd_lcm..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2 - 1
   y = a + 1

   @test gcd(x, y) == F(1)
   @test gcd(F(0), F(0)) == F(0)
end

@testset "n_transExt.Polynomials..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])
   R, (x, ) = PolynomialRing(F, ["x"])

   f = (1 + 2x)*a^2 + 3x*b*c + (x + 1)

   g = f^5

   @test g == (2*a^10-2*b^5*c^5+1)*x^5+(a^10+1)

end
