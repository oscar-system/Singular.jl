@testset "n_Q.constructors" begin

   @test QQ isa Nemo.Field

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

@testset "n_Q.printing" begin
   @test string(QQ(123)) == "123"
   @test sprint(show, "text/plain", (QQ(123))) == "123"
end

@testset "n_Q.manipulation" begin
   @test isone(one(QQ))
   @test iszero(zero(QQ))
   @test is_unit(QQ(1)) && is_unit(QQ(2))
   @test !is_unit(QQ(0))

   @test characteristic(QQ) == 0

   @test numerator(QQ(2)) == 2
   @test denominator(QQ(2)) == 1
   @test denominator(QQ()) == 1

   @test abs(QQ(-2)) == 2
   @test abs(QQ(2)) == 2
   @test abs(QQ()) == 0

   @test deepcopy(QQ(2)) == QQ(2)
end

@testset "n_Q.unary_ops" begin
   @test -QQ(123) == -123
   @test -QQ() == 0
end

@testset "n_Q.binary_ops" begin
   a = QQ(2)
   b = QQ(3)

   @test a + b == 5
   @test a - b == -1
   @test a*b == 6
end

@testset "n_Q.comparison" begin
   @test QQ(2) < QQ(3)
   @test QQ(2) <= QQ(3)
   @test QQ(3) > QQ(2)
   @test QQ(3) >= QQ(3)
   @test QQ(2) == QQ(2)
   @test isequal(QQ(2), QQ(2))
end

@testset "n_Q.ad_hoc_comparison" begin
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

@testset "n_Q.powering" begin
   @test QQ(2)^10 == 1024
end

@testset "n_Q.exact_division" begin
   @test_throws ErrorException inv(zero(QQ))
   @test_throws ErrorException divexact(one(QQ), zero(QQ))
   @test_throws ErrorException divexact(zero(QQ), zero(QQ))

   @test divexact(QQ(12), QQ(2)) == 6
end

@testset "n_Q.gcd_lcm" begin
   @test gcd(QQ(6), QQ(12)) == 1
   @test gcd(QQ(-6), QQ(12)) == 1
   @test gcd(QQ(6), QQ(-12)) == 1
   @test gcd(QQ(-6), QQ(-12)) == 1

   @test gcd(QQ(0), QQ(0)) == 0

   @test lcm(QQ(4), QQ(6)) == 24
end

@testset "n_Q.reconstruct" begin
   @test reconstruct(ZZ(11), ZZ(15)) == divexact(QQ(-1), QQ(4))
   @test reconstruct(ZZ(11), 15) == divexact(QQ(-1), QQ(4))
   @test reconstruct(11, ZZ(15)) == divexact(QQ(-1), QQ(4))
end

@testset "n_Q.Polynomials" begin
   R, x = Nemo.polynomial_ring(QQ, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == 9*x^4+12*x^3+10*x^2+4*x+1
end

@testset "n_Q.conversions" begin
   @test Singular.QQ(1//2) == Singular.QQ(1) // Singular.QQ(2)
   @test AbstractAlgebra.QQ(Singular.QQ(1) // Singular.QQ(2)) == 1//2

   s = Singular.QQ(1//2)
   @test 1//2 == Rational{BigInt}(s) isa Rational{BigInt}
   @test 1//2 == Rational(s) isa Rational{BigInt}
   @test 1//2 == convert(Rational, s) isa Rational{BigInt}
   @test 1//2 == Rational{Int}(s) isa Rational{Int}
   @test 1//2 == convert(Rational{Int}, s) isa Rational{Int}

   @test 1//2 == Nemo.QQ(s) isa Nemo.QQFieldElem
   @test 1//2 == Nemo.QQFieldElem(s) isa Nemo.QQFieldElem
   @test 1//2 == convert(Nemo.QQFieldElem, s) isa Nemo.QQFieldElem
end

@testset "n_Q.rand" begin
   for x in (rand(QQ, 1:9), rand(rng, QQ, 1:9),
             rand(make(QQ, 1:9)))
      @test x isa n_Q
      @test numerator(x) <= 9
      @test denominator(x) <= 9
   end
   Random.seed!(rng, 0)
   M = rand(rng, make(QQ, 1:9), 2, 3)
   @test M isa Matrix{n_Q}
   foreach(M) do x
      @test numerator(x) <= 9
      @test denominator(x) <= 9
   end
   Random.seed!(rng, 0)
   @test M == rand(rng, make(QQ, 1:9), 2, 3)
end

@testset "n_Q.immediate_boundaries" begin
   # Singular stores an integer of up to 29 bits in the number word itself
   for v in (0, 1, -1, 2^28 - 1, 2^28, -2^28, -2^28 - 1, 2^40,
             typemax(Int), typemin(Int), big(2)^100, -big(2)^100)
      a = QQ(v)

      @test a == v && a == Nemo.ZZRingElem(v)
      @test Rational{BigInt}(a) == v
      @test Nemo.QQFieldElem(a) == v
      @test numerator(a) == v && isone(denominator(a))
      @test a + QQ(1) == big(v) + 1
      @test a * QQ(2) == 2 * big(v)

      @test iszero(a) == iszero(v)
      @test isone(a) == isone(v)
      @test is_unit(a) == !iszero(v)
      @test Singular.is_negative(a) == (v < 0)
      @test abs(a) == abs(big(v))
      @test hash(a) == hash(QQ(v))

      # from n_Z, where iszero and isone rely on getting an immediate back
      b = QQ(ZZ(v))
      @test b == a
      @test iszero(b) == iszero(v)
      @test isone(b) == isone(v)
   end
end

@testset "n_Q.nemo_conversions" begin
   for (n, d) in ((3, 7), (-3, 7), (1, 1), (2^30, 7), (7, 2^30),
                  (big(2)^100, 7), (7, big(2)^100), (-big(2)^100, big(3)^50))
      q = Nemo.QQFieldElem(Nemo.ZZRingElem(n), Nemo.ZZRingElem(d))
      a = QQ(q)

      @test a == q
      @test Nemo.QQFieldElem(a) == q
      @test Rational{BigInt}(a) == Rational{BigInt}(big(n), big(d))
      @test numerator(a) == numerator(q)
      @test denominator(a) == denominator(q)
      @test QQ(Rational{BigInt}(q)) == a
      @test hash(a) == hash(QQ(q))
   end
end

@testset "n_Q.fraction_construction" begin
   # QQ(::QQFieldElem) and QQ(::Rational) build the snumber without a gcd,
   # relying on both inputs being in lowest terms with a positive denominator
   for (n, d) in ((3, 7), (-3, 7), (2^28 - 1, 2^28), (-2^28, 2^28 - 1),
                  (typemax(Int), 2), (typemin(Int) + 1, 3),
                  (big(2)^100 + 1, 7), (7, big(3)^50),
                  (-big(2)^100 - 1, big(3)^50))
      ref = QQ(n) // QQ(d)
      for a in (QQ(Nemo.QQ(n, d)), QQ(Rational{BigInt}(n, d)))
         @test a == ref
         @test numerator(a) == n && denominator(a) == d
         @test a * QQ(d) == n
         @test hash(a) == hash(ref)
         @test Singular.libSingular.rational_state(a.ptr) ==
               Singular.libSingular.RAT_NORMALIZED
      end
   end

   @test QQ(3//7) == QQ(3)//QQ(7)
   @test QQ(-3//7) == QQ(-3)//QQ(7)
   @test QQ(5//1) == 5 && QQ(Nemo.QQ(5)) == 5
end

@testset "n_Q.nemo_comparison" begin
   @test QQ(5) == Nemo.ZZ(5)
   @test Nemo.ZZ(5) == QQ(5)
   @test QQ(5) == Nemo.QQ(5)
   @test Nemo.QQ(5) == QQ(5)
   @test QQ(3)//QQ(7) == Nemo.QQ(3, 7)
   @test QQ(3)//QQ(7) != Nemo.QQ(3, 8)
   @test QQ(3)//QQ(7) < Nemo.QQ(1, 2)
   @test Nemo.QQ(1, 2) > QQ(3)//QQ(7)
   @test QQ(3)//QQ(7) < Nemo.ZZ(1)
   @test QQ(3)//QQ(7) == 3//7
   @test 3//7 == QQ(3)//QQ(7)

   # a proper fraction is never an integer, however it is stored
   @test QQ(3)//QQ(7) != 1
   @test QQ(3)//QQ(7) != 0
   @test 0 < QQ(3)//QQ(7) < 1

   @test cmp(QQ(3)//QQ(7), 1) == -1 && cmp(1, QQ(3)//QQ(7)) == 1
   @test cmp(QQ(3)//QQ(7), 3//7) == 0
   @test cmp(QQ(5), Nemo.ZZ(5)) == 0
   @test cmp(Nemo.QQ(1, 2), QQ(3)//QQ(7)) == 1
end

@testset "n_Q.unnormalized_coefficients" begin
   # polynomial arithmetic leaves coefficients unnormalized, i.e. with a
   # numerator and denominator that still share a factor
   R, (x,) = polynomial_ring(QQ, ["x"])
   r = QQ(3)//QQ(7)
   expected = [Nemo.QQ(27, 343), Nemo.QQ(81, 343),
               Nemo.QQ(81, 343), Nemo.QQ(27, 343)]

   for (c, e) in zip(coefficients((r*x + r)^3), expected)
      @test c == e
      @test numerator(c) == numerator(e)
      @test denominator(c) == denominator(e)
      # normalizing replaces the number in place, so repeat every conversion
      @test Nemo.QQFieldElem(c) == e
      @test Nemo.QQFieldElem(c) == e
      @test Rational{BigInt}(c) == Rational{BigInt}(e)
      @test hash(c) == hash(QQ(e))
   end
end
