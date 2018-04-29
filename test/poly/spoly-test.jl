function test_spoly_constructors()
   print("spoly.constructors...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   @test elem_type(R) == spoly{n_Z}
   @test elem_type(PolyRing{n_Z}) == spoly{n_Z}
   @test parent_type(spoly{n_Z}) == PolyRing{n_Z}
   @test base_ring(R) == ZZ

   typeof(R) <: Nemo.Ring

   a = R()

   @test base_ring(a) == ZZ
   @test parent(a) == R

   @test isa(a, spoly)

   b = R(123)

   @test isa(b, spoly)

   c = R(BigInt(123))

   @test isa(c, spoly)

   d = R(c)

   @test isa(d, spoly)

   f = R(Nemo.ZZ(123))

   @test isa(f, spoly)

   g = R(ZZ(123))

   @test isa(g, spoly)

   S, (y, ) = PolynomialRing(QQ, ["y", ])

   h = S(ZZ(123))

   @test isa(h, spoly)

   T, (z, ) = PolynomialRing(Nemo.ZZ, ["z", ])

   k = T(123)

   @test isa(k, spoly)

   S = @PolynomialRing(ZZ, "x", 50)

   @test isa(x17, spoly)

   T = @PolynomialRing(ZZ, "y", 50, :lex)

   @test isa(y7, spoly)

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   @test isgen(x)
   @test isgen(y)
   @test !isgen(R(1))
   @test !isgen(R(0))
   @test !isgen(2x)
   @test !isgen(x + y)
   @test !isgen(x^2)

   println("PASS")
end

function test_spoly_printing()
   print("spoly.printing...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   @test string(3x^2 + 2x + 1) == "3*x^2+2*x+1"

   println("PASS")
end

function test_spoly_manipulation()
   print("spoly.manipulation...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   @test isone(one(R))
   @test iszero(zero(R))
   @test isunit(R(1)) && isunit(R(-1))
   @test !isunit(R(2)) && !isunit(R(0)) && !isunit(x) 
   @test isgen(x)
   @test !isgen(R(1)) && !isgen(x + 1)
   @test isconstant(R(0)) && isconstant(R(1))
   @test !isconstant(x) && !isconstant(x + 1)

   @test length(x^2 + 2x + 1) == 3
   @test degree(x^2 + 2x + 1) == 2

   @test exponent(x^5 + 3x + 2, 2) == [5]

   A = [0]

   exponent!(A, x^5 + 3x + 2, 2)

   @test A == [5]

   @test lead_exponent(x^3 + 2x + 1) == [3]

   @test deepcopy(x + 2) == x + 2

   @test characteristic(R) == 0

   @test ngens(R) == 1

   pol = x^5 + 3x + 2
   c_iter = coeffs_expos(pol)
   
   i = length(pol) - 1
   for (c, ex) in c_iter
      @test coeff(pol, i) == c
      @test exponent(pol, i) == ex
      i -= 1
   end

   R, (x, ) = PolynomialRing(ResidueRing(ZZ, 6), ["x", ])

   @test characteristic(R) == 6

   println("PASS")
end

function test_spoly_unary_ops()
   print("spoly.unary_ops...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1

   @test -a == -x^2 - 3x - 1

   println("PASS")
end

function test_spoly_binary_ops()
   print("spoly.binary_ops...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1
   b = 2x + 4

   @test a + b == x^2+5*x+5
   @test a - b == x^2+x-3
   @test a*b == 2*x^3+10*x^2+14*x+4

   println("PASS")
end

function test_spoly_comparison()
   print("spoly.comparison...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1

   @test a == deepcopy(a)
   @test a != x

   println("PASS")
end

function test_spoly_powering()
   print("spoly.powering...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1

   @test a^0 == 1
   @test a^1 == x^2 + 3x + 1
   @test a^3 == x^6+9*x^5+30*x^4+45*x^3+30*x^2+9*x+1

   println("PASS")
end

function test_spoly_exact_division()
   print("spoly.exact_division...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1
   b = 2x + 4

   @test divexact(a*b, a) == b

   println("PASS")
end

function test_spoly_gcd_lcm()
   print("spoly.gcd_lcm...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1
   b = 2x + 4
   c = 2x^2 + 1

   @test gcd(a*c, b*c) == c

   # segfault in p_Divide
   # @test lcm(a, b) == a*b

   # seems to fail
   # @test primpart(2*a) == a

   @test content(2*a) == 2

   println("PASS")
end

function test_spoly_extended_gcd()
   print("spoly.extended_gcd...")

   R, (x, ) = PolynomialRing(QQ, ["x", ])

   a = x^2 + 3x + 1
   b = 2x + 4

   g, s, t = gcdx(a, b)

   @test s*a + t*b == g

   println("PASS")
end

function test_spoly_Polynomials()
   print("spoly.Polynomials...")

   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   S, y = Nemo.PolynomialRing(R, "y")
   
   f = (1 + 2x + 3x^2)*y + (2x + 3)

   g = f^2

   @test g == (9*x^4+12*x^3+10*x^2+4*x+1)*y^2+(12*x^3+26*x^2+16*x+6)*y+(4*x^2+12*x+9)

   println("PASS")
end

function test_spoly()
   test_spoly_constructors()
   test_spoly_printing()
   test_spoly_manipulation()
   test_spoly_unary_ops()
   test_spoly_binary_ops()
   test_spoly_comparison()
   test_spoly_powering()
   test_spoly_exact_division()
   test_spoly_gcd_lcm()
   test_spoly_extended_gcd()
   test_spoly_Polynomials()

   println("")
end

