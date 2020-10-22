@testset "n_Q.constructors..." begin
   @test elem_type(QQ) == n_Q
   @test elem_type(Rationals) == n_Q
   @test parent_type(n_Q) == Rationals
   @test base_ring(QQ) == ZZ

   typeof(QQ) <: Nemo.Field

   a = QQ()

   @test base_ring(a) == ZZ
   @test parent(a) == QQ

   @test isa(a, n_Q)

   b = QQ(123)

   @test isa(b, n_Q)

   c = QQ(BigInt(123))

   @test isa(c, n_Q)

   d = QQ(c)

   @test isa(d, n_Q)

   f = QQ(Nemo.ZZ(123))

   @test isa(f, n_Q)

   g = QQ(ZZ(123))

   @test isa(g, n_Q)

   h = QQ(1, 2)

   @test isa(h, n_Q)

   h2 = QQ(Nemo.ZZ(2), 3)

   @test isa(h2, n_Q)
   @test h2 == 2//3

   h3 = QQ(ZZ(2), ZZ(3))

   @test isa(h3, n_Q)
   @test h3 == 2//3

   i = QQ(1//2)

   @test isa(i, n_Q)
   @test i == QQ(1) // QQ(2)

   j = QQ(BigInt(1)//BigInt(2))

   @test isa(j, n_Q)
   @test j == QQ(1) // QQ(2)

   k = QQ(Nemo.QQ(2, 3))

   @test isa(k, n_Q)
   @test k == QQ(2) // QQ(3)
end

@testset "n_Q.printing..." begin
   @test string(QQ(123)) == "123"
end

@testset "n_Q.manipulation..." begin
   @test isone(one(QQ))
   @test iszero(zero(QQ))
   @test isunit(QQ(1)) && isunit(QQ(2))
   @test !isunit(QQ(0))

   @test numerator(QQ(2)) == 2
   @test denominator(QQ(2)) == 1
   @test denominator(QQ()) == 1

   @test abs(QQ(-2)) == 2
   @test abs(QQ(2)) == 2
   @test abs(QQ()) == 0

   @test deepcopy(QQ(2)) == QQ(2)
end

@testset "n_Q.unary_ops..." begin
   @test -QQ(123) == -123
   @test -QQ() == 0
end

@testset "n_Q.binary_ops..." begin
   a = QQ(2)
   b = QQ(3)

   @test a + b == 5
   @test a - b == -1
   @test a*b == 6
end

@testset "n_Q.adhoc_binary..." begin
   @test QQ(2) + 3 == 5
   @test 2 + QQ(3) == 5
   @test QQ(2) - 3 == -1
   @test 2 - QQ(3) == -1
   @test 2*QQ(3) == 6
   @test QQ(2)*3 == 6
end

@testset "n_Q.comparison..." begin
   @test QQ(2) < QQ(3)
   @test QQ(2) <= QQ(3)
   @test QQ(3) > QQ(2)
   @test QQ(3) >= QQ(3)
   @test QQ(2) == QQ(2)
   @test isequal(QQ(2), QQ(2))
end

@testset "n_Q.ad_hoc_comparison..." begin
   @test QQ(2) < 3
   @test QQ(2) <= 3
   @test 2 < QQ(3)
   @test 2 <= QQ(3)
   @test QQ(3) > 2
   @test QQ(3) >= 3
   @test 3 > QQ(2)
   @test 3 >= QQ(3)
   @test QQ(2) == 2
   @test 2 == QQ(2)
   @test isequal(QQ(2), 2)
   @test isequal(2, QQ(2))
end

@testset "n_Q.powering..." begin
   @test QQ(2)^10 == 1024
end

@testset "n_Q.exact_division..." begin
   @test divexact(QQ(12), QQ(2)) == 6
end

@testset "n_Q.gcd_lcm..." begin
   @test gcd(QQ(6), QQ(12)) == 1
   @test gcd(QQ(-6), QQ(12)) == 1
   @test gcd(QQ(6), QQ(-12)) == 1
   @test gcd(QQ(-6), QQ(-12)) == 1

   @test gcd(QQ(0), QQ(0)) == 0

   @test lcm(QQ(4), QQ(6)) == 24
end

@testset "n_Q.reconstruct..." begin
   @test reconstruct(ZZ(11), ZZ(15)) == divexact(QQ(-1), QQ(4))
   @test reconstruct(ZZ(11), 15) == divexact(QQ(-1), QQ(4))
   @test reconstruct(11, ZZ(15)) == divexact(QQ(-1), QQ(4))
end

@testset "n_Q.Polynomials..." begin
   R, x = Nemo.PolynomialRing(QQ, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == 9*x^4+12*x^3+10*x^2+4*x+1
end

@testset "n_Q.conversions..." begin
   @test Singular.QQ(1//2) == Singular.QQ(1) // Singular.QQ(2)
   @test AbstractAlgebra.QQ(Singular.QQ(1) // Singular.QQ(2)) == 1//2

   s = Singular.QQ(1//2)
   @test 1//2 == Rational{BigInt}(s) isa Rational{BigInt}
   @test 1//2 == Rational(s) isa Rational{BigInt}
   @test 1//2 == convert(Rational, s) isa Rational{BigInt}
   @test 1//2 == Rational{Int}(s) isa Rational{Int}
   @test 1//2 == convert(Rational{Int}, s) isa Rational{Int}

   @test 1//2 == Nemo.QQ(s) isa Nemo.fmpq
   @test 1//2 == Nemo.fmpq(s) isa Nemo.fmpq
   @test 1//2 == convert(Nemo.fmpq, s) isa Nemo.fmpq
end
