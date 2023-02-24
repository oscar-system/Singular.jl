@testset "n_Zn.constructors" begin
   F = residue_ring(ZZ, 8)
   F1 = residue_ring(ZZ, 8)
   F2 = residue_ring(ZZ, 8, cached = false)

   @test F isa Singular.Ring
   @test F1 isa Singular.Ring
   @test F2 isa Singular.Ring
   @test F == F1
   @test F != F2
   @test F1 != F2

   F = residue_ring(ZZ, BigInt(10)^50)
   F1 = residue_ring(ZZ, BigInt(10)^50)
   F2 = residue_ring(ZZ, BigInt(10)^50, cached = false)

   @test F isa Singular.Ring
   @test F1 isa Singular.Ring
   @test F2 isa Singular.Ring
   @test F == F1
   @test F != F2
   @test F1 != F2

   @test_throws DomainError residue_ring(ZZ, -rand(1:99))

end

@testset "n_Zn.printing" begin
   R = residue_ring(ZZ, 5)

   @test string(R(3)) == "3"

   @test sprint(show, "text/plain", (R(3))) == "3"
end

@testset "n_Zn.manipulation" begin
   R = residue_ring(ZZ, 6)

   @test isone(one(R))
   @test iszero(zero(R))
   @test is_unit(R(1)) && is_unit(R(5))
   @test !is_unit(R(0))  && !is_unit(R(2))

   @test characteristic(R) == 6

   @test deepcopy(R(2)) == R(2)

   NR = Nemo.residue_ring(Nemo.ZZ, 6)
   @test R(2) == R(NR(2))
   @test NR(2) == NR(R(2))

   NR = Nemo.residue_ring(Nemo.ZZ, Nemo.ZZ(6))
   @test R(2) == R(NR(2))
   @test NR(2) == NR(R(2))
end

@testset "n_Zn.unary_ops" begin
   R = residue_ring(ZZ, 5)

   @test -R(3) == R(2)
   @test -R() == R()

   R = residue_ring(ZZ, BigInt(10)^40)
   @test iszero(R(10)^40)

end

@testset "n_Zn.binary_ops" begin
   R = residue_ring(ZZ, 5)

   a = R(2)
   b = R(3)

   @test a + b == R(0)
   @test a - b == R(4)
   @test a*b == R(1)
end

@testset "n_Zn.comparison" begin
   R = residue_ring(ZZ, 5)

   @test R(2) == R(2)
   @test isequal(R(2), R(2))
end

@testset "n_Zn.ad_hoc_comparison" begin
   R = residue_ring(ZZ, 5)

   @test R(2) == 2
   @test 2 == R(2)
   @test isequal(R(2), 2)
   @test isequal(2, R(2))

   @test R(2) == ZZ(2)
   @test ZZ(2) == R(2)
end

@testset "n_Zn.powering" begin
   R = residue_ring(ZZ, 5)

   @test R(2)^10 == R(4)
   @test_throws DomainError R(2)^-rand(1:99)
end

@testset "n_Zn.exact_division" begin
   R = residue_ring(ZZ, 5)

   @test_throws ErrorException inv(zero(R))
   @test_throws ErrorException divexact(one(R), zero(R))
   @test_throws ErrorException divexact(zero(R), zero(R))

   @test divexact(R(2), R(3)) == R(4)
end

@testset "n_Zn.gcd_lcm" begin
   R = residue_ring(ZZ, 5)

   @test gcd(R(2), R(3)) == R(1)
   @test gcd(R(0), R(0)) == R(0)
end

@testset "n_Zn.extended_gcd" begin
   R = residue_ring(ZZ, 6)

   g, s, t = gcdx(R(2), R(4))

   @test g == R(2)*s + R(4)*t

   g, s, t = gcdx(R(1), R(5))

   @test g == R(1)*s + R(5)*t
end

@testset "n_Zn.Polynomials" begin
   R = residue_ring(ZZ, 5)
   S, x = Nemo.polynomial_ring(R, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == -x^4+2*x^3-x+1
end
