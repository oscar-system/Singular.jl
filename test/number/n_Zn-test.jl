function test_n_Zn_constructors()
   print("n_Zn.constructors...")

   R = ResidueRing(ZZ, 7)

   @test elem_type(R) == n_Zn
   @test elem_type(N_ZnRing) == n_Zn
   @test parent_type(n_Zn) == N_ZnRing
   @test base_ring(R) == ZZ

   typeof(R) <: Nemo.Ring

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

   println("PASS")
end

function test_n_Zn_printing()
   print("n_Zn.printing...")

   R = ResidueRing(ZZ, 5)

   @test string(R(3)) == "3"

   println("PASS")
end

function test_n_Zn_manipulation()
   print("n_Zn.manipulation...")

   R = ResidueRing(ZZ, 6)

   @test isone(one(R))
   @test iszero(zero(R))
   @test isunit(R(1)) && isunit(R(5))
   @test !isunit(R(0))  && !isunit(R(2))

   @test characteristic(R) == 6

   @test deepcopy(R(2)) == R(2)

   println("PASS")
end

function test_n_Zn_unary_ops()
   print("n_Zn.unary_ops...")

   R = ResidueRing(ZZ, 5)

   @test -R(3) == R(2)
   @test -R() == R()

   println("PASS")
end

function test_n_Zn_binary_ops()
   print("n_Zn.binary_ops...")

   R = ResidueRing(ZZ, 5)

   a = R(2)
   b = R(3)

   @test a + b == R(0)
   @test a - b == R(4)
   @test a*b == R(1)

   println("PASS")
end

function test_n_Zn_adhoc_binary()
   print("n_Zn.adhoc_binary...")

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

   println("PASS")
end

function test_n_Zn_comparison()
   print("n_Zn.comparison...")

   R = ResidueRing(ZZ, 5)

   @test R(2) == R(2)
   @test isequal(R(2), R(2))

   println("PASS")
end

function test_n_Zn_adhoc_comparison()
   print("n_Zn.ad_hoc_comparison...")

   R = ResidueRing(ZZ, 5)

   @test R(2) == 2
   @test 2 == R(2)
   @test isequal(R(2), 2)
   @test isequal(2, R(2))

   @test R(2) == ZZ(2)
   @test ZZ(2) == R(2)

   println("PASS")
end

function test_n_Zn_powering()
   print("n_Zn.powering...")

   R = ResidueRing(ZZ, 5)

   @test R(2)^10 == R(4)

   println("PASS")
end

function test_n_Zn_exact_division()
   print("n_Zn.exact_division...")

   R = ResidueRing(ZZ, 5)

   @test divexact(R(2), R(3)) == R(4)

   println("PASS")
end

function test_n_Zn_gcd_lcm()
   print("n_Zn.gcd_lcm...")

   R = ResidueRing(ZZ, 5)

   @test gcd(R(2), R(3)) == R(1)
   @test gcd(R(0), R(0)) == R(0)

   println("PASS")
end

function test_n_Zn_extended_gcd()
   print("n_Zn.extended_gcd...")

   R = ResidueRing(ZZ, 6)

   g, s, t = gcdx(R(2), R(4))

   @test g == R(2)*s + R(4)*t

   g, s, t = gcdx(R(1), R(5))

   @test g == R(1)*s + R(5)*t

   println("PASS")
end

function test_n_Zn_Polynomials()
   print("n_Zn.Polynomials...")

   R = ResidueRing(ZZ, 5)
   S, x = Nemo.PolynomialRing(R, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == -x^4+2*x^3-x+1

   println("PASS")
end

function test_n_Zn()
   test_n_Zn_constructors()
   test_n_Zn_printing()
   test_n_Zn_manipulation()
   test_n_Zn_unary_ops()
   test_n_Zn_binary_ops()
   test_n_Zn_adhoc_binary()
   test_n_Zn_comparison()
   test_n_Zn_adhoc_comparison()
   test_n_Zn_powering()
   test_n_Zn_exact_division()
   test_n_Zn_gcd_lcm()
   test_n_Zn_extended_gcd()
   test_n_Zn_Polynomials()

   println("")
end

