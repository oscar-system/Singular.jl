@testset "n_FF.constructors..." begin
   F, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])
   
   @test elem_type(F) == n_FF
   @test elem_type(N_FField) == n_FF
   @test parent_type(n_FF) == N_FField
   @test base_ring(F) == QQ

   typeof(F) <: Singular.Field

   a = F()

   @test base_ring(a) == Union{}
   @test parent(a) == F

   @test isa(a, n_FF)

   b = F(3)

   @test isa(b, n_FF)

   c = F(BigInt(3))

   @test isa(c, n_FF)

   # segfault
   # d = R(ZZ(3))

   # @test isa(d, n_FF);

   f = F(c)

   @test isa(f, n_FF)

   f = F(Nemo.ZZ(123))

   @test isa(f, n_FF)
end

@testset "n_FF.printing..." begin
   F, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])

   @test string(3*a*b + 2*c) == "(3*a*b+2*c)"
end

@testset "n_FF.manipulation..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   @test isone(one(F))
   @test iszero(zero(F))
   @test isunit(F(1)) && isunit(F(2))
   @test !isunit(F(0)) 

   @test characteristic(F) == 5
   @test transcendence_degree(F) == 3

   @test deepcopy(F(2)) == F(2)
end

@testset "n_FF.unary_ops..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = 3*a + b

   @test -F(3) == F(2)
   @test -F() == F()
   @test -x == 2*a - b
end

@testset "n_FF.binary_ops..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = 2*a + 2*c
   y = 3*b + 2*c

   @test x + y == 2*a + 3*b - c
   @test x - y == 2*a + 2*b
   @test x*y == a*b - a*c + b*c - c^2
end

@testset "n_FF.comparison..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2*b*c + b^2*c + 5*a

   @test x == deepcopy(x)
   @test isequal(x, x)
end

@testset "n_FF.powering..." begin
   F, (a, b, c) = FunctionField(Fp(3), ["a", "b", "c"])

   @test (a*b*c + 1)^3 == a^3*b^3*c^3 + 1
end

@testset "n_FF.exact_division..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2 - 1
   y = a + 1

   @test inv(x) == 1 // (a^2 - 1)
   @test divexact(x, y) == a-1
end

@testset "n_FF.gcd_lcm..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])

   x = a^2 - 1
   y = a + 1

   @test gcd(x, y) == F(1)
   @test gcd(F(0), F(0)) == F(0)
end

@testset "n_FF.Polynomials..." begin
   F, (a, b, c) = FunctionField(Fp(5), ["a", "b", "c"])
   R, (x, ) = PolynomialRing(F, ["x"])

   f = (1 + 2x)*a^2 + 3x*b*c + (x + 1)

   g = f^5

   @test g == (2*a^10-2*b^5*c^5+1)*x^5+(a^10+1)

end

