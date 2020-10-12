@testset "n_Z.constructors..." begin
   @test elem_type(ZZ) == n_Z
   @test elem_type(Integers) == n_Z
   @test parent_type(n_Z) == Integers
   @test base_ring(ZZ) == Union{}

   typeof(ZZ) <: Nemo.Ring

   a = ZZ()

   @test base_ring(a) == Union{}
   @test parent(a) == ZZ

   @test isa(a, n_Z)

   b = ZZ(123)

   @test isa(b, n_Z)

   c = ZZ(BigInt(123))

   @test isa(c, n_Z)

   d = ZZ(c)

   @test isa(d, n_Z)

   f = ZZ(Nemo.ZZ(123))

   # should not crash
   for i=1:100
      f = ZZ(Nemo.ZZ(3)^(2*i)-1)
   end

   @test isa(f, n_Z)
end

@testset "n_Z.printing..." begin
   @test string(ZZ(123)) == "123"
end

@testset "n_Z.manipulation..." begin
   @test isone(one(ZZ))
   @test iszero(zero(ZZ))
   @test isunit(ZZ(1)) && isunit(ZZ(-1))
   @test !isunit(ZZ(2)) && !isunit(ZZ(0))

   @test numerator(ZZ(2)) == 2
   @test denominator(ZZ(2)) == 1
   @test denominator(ZZ()) == 1

   @test abs(ZZ(-2)) == 2
   @test abs(ZZ(2)) == 2
   @test abs(ZZ()) == 0

   @test deepcopy(ZZ(2)) == ZZ(2)
end

@testset "n_Z.unary_ops..." begin
   @test -ZZ(123) == -123
   @test -ZZ() == 0
end

@testset "n_Z.binary_ops..." begin
   a = ZZ(2)
   b = ZZ(3)

   @test a + b == 5
   @test a - b == -1
   @test a*b == 6
end

@testset "n_Z.adhoc_binary..." begin
   @test ZZ(2) + 3 == 5
   @test 2 + ZZ(3) == 5
   @test ZZ(2) - 3 == -1
   @test 2 - ZZ(3) == -1
   @test 2*ZZ(3) == 6
   @test ZZ(2)*3 == 6
end

@testset "n_Z.comparison..." begin
   @test ZZ(2) < ZZ(3)
   @test ZZ(2) <= ZZ(3)
   @test ZZ(3) > ZZ(2)
   @test ZZ(3) >= ZZ(3)
   @test ZZ(2) == ZZ(2)
   @test isequal(ZZ(2), ZZ(2))
end

@testset "n_Z.ad_hoc_comparison..." begin
   @test ZZ(2) < 3
   @test ZZ(2) <= 3
   @test 2 < ZZ(3)
   @test 2 <= ZZ(3)
   @test ZZ(3) > 2
   @test ZZ(3) >= 3
   @test 3 > ZZ(2)
   @test 3 >= ZZ(3)
   @test ZZ(2) == 2
   @test 2 == ZZ(2)
   @test isequal(ZZ(2), 2)
   @test isequal(2, ZZ(2))

   # fmpz
   @test ZZ(3) == Nemo.fmpz(3)
   @test Nemo.fmpz(3) == ZZ(3)
   @test ZZ(3)^42 == Nemo.fmpz(3)^42
   @test Nemo.fmpz(3)^42 == ZZ(3)^42
end

@testset "n_Z.powering..." begin
   @test ZZ(2)^10 == 1024
end

@testset "n_Z.exact_division..." begin
   @test divexact(ZZ(12), ZZ(2)) == 6
end

@testset "n_Z.euclidean_division..." begin
   @test div(ZZ(7), ZZ(3)) == 2
   @test rem(ZZ(4), ZZ(3)) == 1
   @test mod(ZZ(7), ZZ(3)) == 1
   @test rem(ZZ(-2), ZZ(3)) == -2
   @test rem(ZZ(2), ZZ(-3)) == 2
   @test rem(ZZ(-2), ZZ(-3)) == -2
   @test mod(ZZ(-2), ZZ(3)) == 1
   @test mod(ZZ(2), ZZ(-3)) == 2
   @test mod(ZZ(-2), ZZ(-3)) == 1
end

@testset "n_Z.gcd_lcm..." begin
   @test gcd(ZZ(6), ZZ(12)) == 6
   @test gcd(ZZ(-6), ZZ(12)) == 6
   @test gcd(ZZ(6), ZZ(-12)) == 6
   @test gcd(ZZ(-6), ZZ(-12)) == 6

   @test lcm(ZZ(4), ZZ(6)) == 12
end

@testset "n_Z.extended_gcd..." begin
   g, s, t = gcdx(ZZ(4), ZZ(6))

   @test s*ZZ(4) + t*ZZ(6) == g
end

@testset "n_Z.chinese_remainder..." begin
#    @test crt(ZZ(2), ZZ(3), ZZ(3), ZZ(7), true) == -4
#    @test crt(ZZ(2), ZZ(3), ZZ(3), ZZ(7), false) == 17
end

@testset "n_Z.Polynomials..." begin
   R, x = Nemo.PolynomialRing(ZZ, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == 9*x^4+12*x^3+10*x^2+4*x+1
end

@testset "n_Z.conversions..." begin
   for n in [-1, 0, 1, -BigInt(2)^65, BigInt(2)^65]
      N = ZZ(n)
      @test convert(BigInt, N) == n
      @test BigInt(N) == n
      for z = (convert(Integer, N), Integer(N))
         @test z isa BigInt
         @test z == n
      end
      @test convert(Int128, N) == n
      @test Int128(N) == n
      @test ZZ(Nemo.fmpz(N)) == N
      @test ZZ(Nemo.ZZ(N)) == N
      @test ZZ(convert(Nemo.fmpz, N)) == N
      if n >= 0
         @test convert(UInt128, N) == n
         @test UInt128(N) == n
      else
         @test_throws InexactError convert(UInt128, N)
      end
      if log2(abs(n)) < 63
         @test convert(Int64, N) == n
         @test Int64(N) == n
         @test convert(Int, N) == n
         @test Int(N) == n
         if n >= 0
            @test convert(UInt64, N) == n
            @test convert(UInt, N) == n
            @test UInt64(N) == n
            @test UInt(N) == n
         end
      else
         @test_throws InexactError convert(Int64, N)
         @test_throws InexactError convert(UInt64, N)
      end

      b = BigInt(N)
      # check no aliasing
      @test b !== BigInt(N)
      b.size = 0
      if n != 0
         @test N != 0
      end
   end
end
