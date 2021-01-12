@testset "n_Zn.constructors..." begin
   R = ResidueRing(ZZ, 7)

   @test elem_type(R) == n_Zn
   @test elem_type(N_ZnRing) == n_Zn
   @test parent_type(n_Zn) == N_ZnRing
   @test base_ring(R) == ZZ

   @test typeof(R) <: Nemo.Ring

   a = R()

   @test base_ring(a) == ZZ
   @test parent(a) == R

   @test isa(a, n_Zn)

   b = R(3)

   @test isa(b, n_Zn)

   c = R(BigInt(3))

   @test isa(c, n_Zn)

   d = R(ZZ(3))

   @test isa(d, n_Zn);

   f = R(c)

   @test isa(f, n_Zn)

   f = R(Nemo.ZZ(123))

   @test isa(f, n_Zn)

   @test_throws DomainError ResidueRing(ZZ, -rand(1:99))
end

@testset "n_Zn.printing..." begin
   R = ResidueRing(ZZ, 5)

   @test string(R(3)) == "3"

   @test sprint(show, "text/plain", (R(3))) == "3"
end

@testset "n_Zn.manipulation..." begin
   R = ResidueRing(ZZ, 6)

   @test isone(one(R))
   @test iszero(zero(R))
   @test isunit(R(1)) && isunit(R(5))
   @test !isunit(R(0))  && !isunit(R(2))

   @test characteristic(R) == 6

   @test deepcopy(R(2)) == R(2)
end

@testset "n_Zn.unary_ops..." begin
   R = ResidueRing(ZZ, 5)

   @test -R(3) == R(2)
   @test -R() == R()
end

@testset "n_Zn.binary_ops..." begin
   R = ResidueRing(ZZ, 5)

   a = R(2)
   b = R(3)

   @test a + b == R(0)
   @test a - b == R(4)
   @test a*b == R(1)
end

@testset "n_Zn.adhoc_binary..." begin
   R = ResidueRing(ZZ, 5)

   @test R(2) + 3 == R(0)
   @test 2 + R(3) == R(0)
   @test R(2) - 3 == R(4)
   @test 2 - R(3) == R(4)
   @test 2*R(3) == R(1)
   @test R(2)*3 == R(1)
   @test R(2) + ZZ(3) == R(0)
   @test ZZ(2) + R(3) == R(0)
   @test R(2) - ZZ(3) == R(4)
   @test ZZ(2) - R(3) == R(4)
   @test ZZ(2)*R(3) == R(1)
   @test R(2)*ZZ(3) == R(1)
end

@testset "n_Zn.comparison..." begin
   R = ResidueRing(ZZ, 5)

   @test R(2) == R(2)
   @test isequal(R(2), R(2))
end

@testset "n_Zn.ad_hoc_comparison..." begin
   R = ResidueRing(ZZ, 5)

   @test R(2) == 2
   @test 2 == R(2)
   @test isequal(R(2), 2)
   @test isequal(2, R(2))

   @test R(2) == ZZ(2)
   @test ZZ(2) == R(2)
end

@testset "n_Zn.powering..." begin
   R = ResidueRing(ZZ, 5)

   @test R(2)^10 == R(4)
   @test_throws DomainError R(2)^-rand(1:99)
end

@testset "n_Zn.exact_division..." begin
   R = ResidueRing(ZZ, 5)

   @test divexact(R(2), R(3)) == R(4)
end

@testset "n_Zn.gcd_lcm..." begin
   R = ResidueRing(ZZ, 5)

   @test gcd(R(2), R(3)) == R(1)
   @test gcd(R(0), R(0)) == R(0)
end

@testset "n_Zn.extended_gcd..." begin
   R = ResidueRing(ZZ, 6)

   g, s, t = gcdx(R(2), R(4))

   @test g == R(2)*s + R(4)*t

   g, s, t = gcdx(R(1), R(5))

   @test g == R(1)*s + R(5)*t
end

@testset "n_Zn.Polynomials..." begin
   R = ResidueRing(ZZ, 5)
   S, x = Nemo.PolynomialRing(R, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == -x^4+2*x^3-x+1
end
