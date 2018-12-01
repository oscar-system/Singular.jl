function test_nemo_fmpq()
   print("Nemo.fmpq...")

   R, (x, y) = PolynomialRing(Nemo.QQ, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   @test isa(coeff(f1, 1), Singular.n_unknown{Nemo.fmpq})

   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2

   @test f1 + Nemo.QQ(2) == Nemo.QQ(2) + f1
   @test f1 - Nemo.QQ(2) == -(Nemo.QQ(2) - f1)
   @test Nemo.QQ(2)*f1 == f1*Nemo.QQ(2)

   @test f1*x == x*f1

   @test deepcopy(f1) == f1

   @test f1*f2 == f2*f1

   @test divexact(f1*f2, f1) == f2

   @test f1^3*f1^3 == f1^6

   @test inv(coeff(f1, 1)) == Nemo.QQ(1, 3)

   @test gcd(coeff(f1, 1), coeff(f1, 2)) == Nemo.QQ(1)

   @test divexact(coeff(f1, 1), coeff(f1, 2)) == Nemo.QQ(3)

   @test coeff(f1, 1) - coeff(f1, 2) == Nemo.QQ(2)

   @test coeff(f1, 1) + coeff(f1, 2) == Nemo.QQ(4)

   println("PASS")
end

function test_nemo_fmpz()
   print("Nemo.fmpz...")

   R, (x, y) = PolynomialRing(Nemo.ZZ, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   @test isa(coeff(f1, 1), Singular.n_unknown{Nemo.fmpz})

   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2

   @test f1 + Nemo.ZZ(2) == Nemo.ZZ(2) + f1
   @test f1 - Nemo.ZZ(2) == -(Nemo.ZZ(2) - f1)
   @test Nemo.ZZ(2)*f1 == f1*Nemo.ZZ(2)

   @test f1*x == x*f1

   @test deepcopy(f1) == f1

   @test f1*f2 == f2*f1

   @test divexact(f1*f2, f1) == f2

   @test f1^3*f1^3 == f1^6

   @test gcd(coeff(f1, 1), coeff(f1, 2)) == Nemo.ZZ(1)

   @test divexact(coeff(f1, 1), coeff(f1, 2)) == Nemo.ZZ(3)

   @test coeff(f1, 1) - coeff(f1, 2) == Nemo.ZZ(2)

   @test coeff(f1, 1) + coeff(f1, 2) == Nemo.ZZ(4)

   println("PASS")
end

function test_nemo_fq_nmod()
   print("Nemo.fq_nmod...")

   F, a = Nemo.FiniteField(7, 2, "a")

   R, (x, y) = PolynomialRing(F, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   @test isa(coeff(f1, 1), Singular.n_unknown{Nemo.fq_nmod})

   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2

   @test f1 + F(2) == F(2) + f1
   @test f1 - F(2) == -(F(2) - f1)
   @test F(2)*f1 == f1*F(2)

   @test f1*x == x*f1

   @test f1 + a == a + f1
   @test f1 - a == -(a - f1)
   @test a*f1 == f1*a

   @test deepcopy(f1) == f1

   @test f1*f2 == f2*f1

   @test divexact(f1*f2, f1) == f2

   @test f1^3*f1^3 == f1^6

   @test inv(coeff(f1, 1)) == F(5)

   @test gcd(coeff(f1, 1), coeff(f1, 2)) == F(1)

   @test divexact(coeff(f1, 1), coeff(f1, 2)) == F(3)

   @test coeff(f1, 1) - coeff(f1, 2) == F(2)

   @test coeff(f1, 1) + coeff(f1, 2) == F(4)

   println("PASS")
end

function test_nemo_fq()
   print("Nemo.fq...")

   F, a = Nemo.FiniteField(Nemo.ZZ(7), 2, "a")

   R, (x, y) = PolynomialRing(F, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   @test isa(coeff(f1, 1), Singular.n_unknown{Nemo.fq})

   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2

   @test f1 + F(2) == F(2) + f1
   @test f1 - F(2) == -(F(2) - f1)
   @test F(2)*f1 == f1*F(2)

   @test f1*x == x*f1

   @test f1 + a == a + f1
   @test f1 - a == -(a - f1)
   @test a*f1 == f1*a

   @test deepcopy(f1) == f1

   @test f1*f2 == f2*f1

   @test divexact(f1*f2, f1) == f2

   @test f1^3*f1^3 == f1^6

   @test inv(coeff(f1, 1)) == F(5)

   @test gcd(coeff(f1, 1), coeff(f1, 2)) == F(1)

   @test divexact(coeff(f1, 1), coeff(f1, 2)) == F(3)

   @test coeff(f1, 1) - coeff(f1, 2) == F(2)

   @test coeff(f1, 1) + coeff(f1, 2) == F(4)

   println("PASS")
end

function test_nemo_nf_elem()
   print("Nemo.nf_elem...")

   U, z = Nemo.PolynomialRing(Nemo.QQ, "z")
   K, a = Nemo.NumberField(z^3 + 3z + 1, "a")

   R, (x, y) = PolynomialRing(K, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   @test isa(coeff(f1, 1), Singular.n_unknown{Nemo.nf_elem})

   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2

   @test f1 + K(2) == K(2) + f1
   @test f1 - K(2) == -(K(2) - f1)
   @test K(2)*f1 == f1*K(2)

   @test f1 + a == a + f1
   @test f1 - a == -(a - f1)
   @test a*f1 == f1*a

   @test f1*x == x*f1

   @test deepcopy(f1) == f1

   @test f1*f2 == f2*f1

   @test divexact(f1*f2, f1) == f2

   @test f1^3*f1^3 == f1^6

   @test inv(coeff(f1, 1)) == K(1)//3

   @test gcd(coeff(f1, 1), coeff(f1, 2)) == K(1)

   @test divexact(coeff(f1, 1), coeff(f1, 2)) == K(3)

   @test coeff(f1, 1) - coeff(f1, 2) == K(2)

   @test coeff(f1, 1) + coeff(f1, 2) == K(4)

   println("PASS")
end

function test_nemo_field()
   print("Nemo.NemoField...")

   U = Nemo.Generic.FractionField(Nemo.ZZ)

   R, (x, y) = PolynomialRing(U, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   @test isa(coeff(f1, 1), Singular.n_unknown{AbstractAlgebra.Generic.Frac{Nemo.fmpz}})

   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2

   @test f1 + U(2) == U(2) + f1
   @test f1 - U(2) == -(U(2) - f1)
   @test U(2)*f1 == f1*U(2)

   @test f1*x == x*f1

   @test deepcopy(f1) == f1

   @test f1*f2 == f2*f1

   @test divexact(f1*f2, f1) == f2

   @test f1^3*f1^3 == f1^6

   @test inv(coeff(f1, 1)) == U(1, 3)

   @test gcd(coeff(f1, 1), coeff(f1, 2)) == U(1)

   @test divexact(coeff(f1, 1), coeff(f1, 2)) == U(3)

   @test coeff(f1, 1) - coeff(f1, 2) == U(2)

   @test coeff(f1, 1) + coeff(f1, 2) == U(4)

   println("PASS")
end

function test_nemo_ring()
   print("Nemo.NemoRing...")

   U, z = Nemo.PolynomialRing(Nemo.ZZ, "z")

   R, (x, y) = PolynomialRing(U, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   @test isa(coeff(f1, 1), Singular.n_unknown{Nemo.fmpz_poly})

   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2

   @test f1 + U(2) == U(2) + f1
   @test f1 - U(2) == -(U(2) - f1)
   @test U(2)*f1 == f1*U(2)

   @test f1*x == x*f1

   @test deepcopy(f1) == f1

   @test f1*f2 == f2*f1

   @test divexact(f1*f2, f1) == f2

   @test f1^3*f1^3 == f1^6

   @test gcd(coeff(f1, 1), coeff(f1, 2)) == U(1)

   @test divexact(coeff(f1, 1), coeff(f1, 2)) == U(3)

   @test coeff(f1, 1) - coeff(f1, 2) == U(2)

   @test coeff(f1, 1) + coeff(f1, 2) == U(4)

   println("PASS")
end

function test_nemo()
#    test_nemo_fmpq()
   test_nemo_fmpz()
#    test_nemo_fq_nmod()
#    test_nemo_fq()
#    test_nemo_nf_elem()
#    test_nemo_field()
#    test_nemo_ring()

   println("")
end
