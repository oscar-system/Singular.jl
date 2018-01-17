function test_n_Zp_constructors()
   print("n_Zp.constructors...")

   R = Fp(7)

   @test elem_type(R) == n_Zp
   @test elem_type(N_ZpField) == n_Zp
   @test parent_type(n_Zp) == N_ZpField
   @test base_ring(R) == Union{}

   typeof(R) <: Nemo.Field

   a = R()

   @test base_ring(a) == Union{}
   @test parent(a) == R

   @test isa(a, n_Zp)

   b = R(3)

   @test isa(b, n_Zp)

   c = R(BigInt(3))

   @test isa(c, n_Zp)

   d = R(ZZ(3))

   @test isa(d, n_Zp);

   f = R(c)

   @test isa(f, n_Zp)

   f = R(Nemo.ZZ(123))

   @test isa(f, n_Zp)

   println("PASS")
end

function test_n_Zp_printing()
   print("n_Zp.printing...")

   @test string(ZZ(123)) == "123"

   println("PASS")
end

function test_n_Zp_manipulation()
   print("n_Zp.manipulation...")

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

function test_n_Zp_unary_ops()
   print("n_Zp.unary_ops...")

   @test -ZZ(123) == -123
   @test -ZZ() == 0

   println("PASS")
end

function test_n_Zp_binary_ops()
   print("n_Zp.binary_ops...")

   a = ZZ(2)
   b = ZZ(3)

   @test a + b == 5
   @test a - b == -1
   @test a*b == 6

   println("PASS")
end

function test_n_Zp_adhoc_binary()
   print("n_Zp.adhoc_binary...")

   @test ZZ(2) + 3 == 5
   @test 2 + ZZ(3) == 5
   @test ZZ(2) - 3 == -1
   @test 2 - ZZ(3) == -1
   @test 2*ZZ(3) == 6
   @test ZZ(2)*3 == 6

   println("PASS")
end

function test_n_Zp_comparison()
   print("n_Zp.comparison...")

   @test ZZ(2) < ZZ(3)
   @test ZZ(2) <= ZZ(3)
   @test ZZ(3) > ZZ(2)
   @test ZZ(3) >= ZZ(3)
   @test ZZ(2) == ZZ(2)
   @test isequal(ZZ(2), ZZ(2))

   println("PASS")
end

function test_n_Zp_adhoc_comparison()
   print("n_Zp.ad_hoc_comparison...")

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

function test_n_Zp_powering()
   print("n_Zp.powering...")

   @test ZZ(2)^10 == 1024

   println("PASS")
end

function test_n_Zp_exact_division()
   print("n_Zp.exact_division...")

   @test divexact(ZZ(12), ZZ(2)) == 6

   println("PASS")
end

function test_n_Zp_euclidean_division()
   print("n_Zp.euclidean_division...")

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

function test_n_Zp_gcd_lcm()
   print("n_Zp.gcd_lcm...")

   @test gcd(ZZ(6), ZZ(12)) == 6
   @test gcd(ZZ(-6), ZZ(12)) == 6
   @test gcd(ZZ(6), ZZ(-12)) == 6
   @test gcd(ZZ(-6), ZZ(-12)) == 6

   @test lcm(ZZ(4), ZZ(6)) == 12

   println("PASS")
end

function test_n_Zp_extended_gcd()
   print("n_Zp.extended_gcd...")

   g, s, t = gcdx(ZZ(4), ZZ(6))

   @test s*ZZ(4) + t*ZZ(6) == g

   println("PASS")
end

function test_n_Zp_chinese_remainder()
   print("n_Zp.chinese_remainder...")

   # @test crt(ZZ(2), ZZ(3), ZZ(3), ZZ(7), true) == -4
   # @test crt(ZZ(2), ZZ(3), ZZ(3), ZZ(7), false) == 17

   println("PASS")
end

function test_n_Zp_Polynomials()
   print("n_Zp.Polynomials...")

   R, x = Nemo.PolynomialRing(ZZ, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == 9*x^4+12*x^3+10*x^2+4*x+1

   println("PASS")
end

function test_n_Zp()
   test_n_Zp_constructors()
   test_n_Zp_printing()
   test_n_Zp_manipulation()
   test_n_Zp_unary_ops()
   test_n_Zp_binary_ops()
   test_n_Zp_adhoc_binary()
   test_n_Zp_comparison()
   test_n_Zp_adhoc_comparison()
   test_n_Zp_powering()
   test_n_Zp_exact_division()
   test_n_Zp_euclidean_division()
   test_n_Zp_gcd_lcm()
   test_n_Zp_extended_gcd()
   test_n_Zp_chinese_remainder()
   test_n_Zp_Polynomials()

   println("")
end

