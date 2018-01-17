function test_n_Z_constructors()
   print("n_Z.constructors...")

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

   @test isa(f, n_Z)

   println("PASS")
end

function test_n_Z_printing()
   print("n_Z.printing...")

   @test string(ZZ(123)) == "123"

   println("PASS")
end

function test_n_Z_manipulation()
   print("n_Z.manipulation...")

   @test isone(one(ZZ))
   @test iszero(zero(ZZ))
   @test isunit(ZZ(1)) && isunit(ZZ(1))
   @test !isunit(ZZ(2)) && !isunit(ZZ(0)) 

   @test numerator(ZZ(2)) == 2
   @test denominator(ZZ(2)) == 1
   @test denominator(ZZ()) == 1

   @test abs(ZZ(-2)) == 2
   @test abs(ZZ(2)) == 2
   @test abs(ZZ()) == 0

   @test deepcopy(ZZ(2)) == ZZ(2)

   println("PASS")
end

function test_n_Z_unary_ops()
   print("n_Z.unary_ops...")

   @test -ZZ(123) == -123
   @test -ZZ() == 0

   println("PASS")
end

function test_n_Z_binary_ops()
   print("n_Z.binary_ops...")

   a = ZZ(2)
   b = ZZ(3)

   @test a + b == 5
   @test a - b == -1
   @test a*b == 6

   println("PASS")
end

function test_n_Z_adhoc_binary()
   print("n_Z.adhoc_binary...")

   @test ZZ(2) + 3 == 5
   @test 2 + ZZ(3) == 5
   @test ZZ(2) - 3 == -1
   @test 2 - ZZ(3) == -1
   @test 2*ZZ(3) == 6
   @test ZZ(2)*3 == 6

   println("PASS")
end

function test_n_Z_comparison()
   print("n_Z.comparison...")

   @test ZZ(2) < ZZ(3)
   @test ZZ(2) <= ZZ(3)
   @test ZZ(3) > ZZ(2)
   @test ZZ(3) >= ZZ(3)
   @test ZZ(2) == ZZ(2)
   @test isequal(ZZ(2), ZZ(2))

   println("PASS")
end

function test_n_Z_adhoc_comparison()
   print("n_Z.ad_hoc_comparison...")

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

   println("PASS")
end

function test_n_Z_powering()
   print("n_Z.powering...")

   @test ZZ(2)^10 == 1024

   println("PASS")
end

function test_n_Z_exact_division()
   print("n_Z.exact_division...")

   @test divexact(ZZ(12), ZZ(2)) == 6

   println("PASS")
end

function test_n_Z_euclidean_division()
   print("n_Z.euclidean_division...")

   @test div(ZZ(7), ZZ(3)) == 2
   @test rem(ZZ(4), ZZ(3)) == 1
   @test mod(ZZ(7), ZZ(3)) == 1
   @test rem(ZZ(-2), ZZ(3)) == -2
   @test rem(ZZ(2), ZZ(-3)) == 2
   @test rem(ZZ(-2), ZZ(-3)) == -2
   @test mod(ZZ(-2), ZZ(3)) == 1
   @test mod(ZZ(2), ZZ(-3)) == 2
   @test mod(ZZ(-2), ZZ(-3)) == 1

   println("PASS")
end

function test_n_Z_gcd_lcm()
   print("n_Z.gcd_lcm...")

   @test gcd(ZZ(6), ZZ(12)) == 6
   @test gcd(ZZ(-6), ZZ(12)) == 6
   @test gcd(ZZ(6), ZZ(-12)) == 6
   @test gcd(ZZ(-6), ZZ(-12)) == 6

   @test lcm(ZZ(4), ZZ(6)) == 12

   println("PASS")
end

function test_n_Z_extended_gcd()
   print("n_Z.extended_gcd...")

   g, s, t = gcdx(ZZ(4), ZZ(6))

   @test s*ZZ(4) + t*ZZ(6) == g

   println("PASS")
end

function test_n_Z_chinese_remainder()
   print("n_Z.chinese_remainder...")

   # @test crt(ZZ(2), ZZ(3), ZZ(3), ZZ(7), true) == -4
   # @test crt(ZZ(2), ZZ(3), ZZ(3), ZZ(7), false) == 17

   println("PASS")
end

function test_n_Z_Polynomials()
   print("n_Z.Polynomials...")

   R, x = Nemo.PolynomialRing(ZZ, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == 9*x^4+12*x^3+10*x^2+4*x+1

   println("PASS")
end

function test_n_Z()
   test_n_Z_constructors()
   test_n_Z_printing()
   test_n_Z_manipulation()
   test_n_Z_unary_ops()
   test_n_Z_binary_ops()
   test_n_Z_adhoc_binary()
   test_n_Z_comparison()
   test_n_Z_adhoc_comparison()
   test_n_Z_powering()
   test_n_Z_exact_division()
   test_n_Z_euclidean_division()
   test_n_Z_gcd_lcm()
   test_n_Z_extended_gcd()
   test_n_Z_chinese_remainder()
   test_n_Z_Polynomials()

   println("")
end

