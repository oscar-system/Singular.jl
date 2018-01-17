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

   R = Fp(5)

   @test string(R(3)) == "-2"

   println("PASS")
end

function test_n_Zp_manipulation()
   print("n_Zp.manipulation...")

   R = Fp(5)

   @test isone(one(R))
   @test iszero(zero(R))
   @test isunit(R(1)) && isunit(R(2))
   @test !isunit(R(0)) 

   @test characteristic(R) == 5

   @test deepcopy(R(2)) == R(2)

   @test Int(R(2)) == 2

   println("PASS")
end

function test_n_Zp_unary_ops()
   print("n_Zp.unary_ops...")

   R = Fp(5)

   @test -R(3) == R(2)
   @test -R() == R()

   println("PASS")
end

function test_n_Zp_binary_ops()
   print("n_Zp.binary_ops...")

   R = Fp(5)

   a = R(2)
   b = R(3)

   @test a + b == R(0)
   @test a - b == R(4)
   @test a*b == R(1)

   println("PASS")
end

function test_n_Zp_adhoc_binary()
   print("n_Zp.adhoc_binary...")

   R = Fp(5)

   @test R(2) + 3 == R(0)
   @test 2 + R(3) == R(0)
   @test R(2) - 3 == R(4)
   @test 2 - R(3) == R(4)
   @test 2*R(3) == R(1)
   @test R(2)*3 == R(1)

   println("PASS")
end

function test_n_Zp_comparison()
   print("n_Zp.comparison...")

   R = Fp(5)

   @test R(2) == R(2)
   @test isequal(R(2), R(2))

   println("PASS")
end

function test_n_Zp_adhoc_comparison()
   print("n_Zp.ad_hoc_comparison...")

   R = Fp(5)

   @test R(2) == 2
   @test 2 == R(2)
   @test isequal(R(2), 2)
   @test isequal(2, R(2))

   println("PASS")
end

function test_n_Zp_powering()
   print("n_Zp.powering...")

   R = Fp(5)

   @test R(2)^10 == R(4)

   println("PASS")
end

function test_n_Zp_exact_division()
   print("n_Zp.exact_division...")

   R = Fp(5)

   @test inv(R(2)) == R(3)
   @test divexact(R(2), R(3)) == R(4)

   println("PASS")
end

function test_n_Zp_gcd_lcm()
   print("n_Zp.gcd_lcm...")

   R = Fp(5)

   @test gcd(R(2), R(3)) == R(1)
   @test gcd(R(0), R(0)) == R(0)

   println("PASS")
end

function test_n_Zp_Polynomials()
   print("n_Zp.Polynomials...")

   R = Fp(5)
   S, x = Nemo.PolynomialRing(R, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == -x^4+2*x^3-x+1

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
   test_n_Zp_gcd_lcm()
   test_n_Zp_Polynomials()

   println("")
end

