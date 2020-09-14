@testset "Nemo.fmpq..." begin
   R, (x, y) = PolynomialRing(Nemo.QQ, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coeffs(f1)]
   @test isa(f1c[1], Singular.n_unknown{Nemo.fmpq})

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

   @test inv(f1c[2]) == Nemo.QQ(1, 3)

   @test gcd(f1c[1], f1c[2]) == Nemo.QQ(1)

   @test divexact(f1c[2], f1c[1]) == Nemo.QQ(3)

   @test f1c[2] - f1c[1] == Nemo.QQ(2)

   @test f1c[1] + f1c[2] == Nemo.QQ(4)

   @test std(Ideal(R, x*y-1, x^2))[1] == Nemo.QQ(1)

   @test length(string((x+y)^2)) > 3

   @test hash((x+y)^2) == hash(x^2+2*x*y+y^2)

   @test deepcopy(f1c[1]) == f1c[1]

   @test canonical_unit(f1c[1]) != 0
end

@testset "Nemo.fmpz..." begin
   R, (x, y) = PolynomialRing(Nemo.ZZ, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coeffs(f1)]

   @test isa(f1c[1], Singular.n_unknown{Nemo.fmpz})

   @test string(first(coeffs(f3))) == "1"

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

   @test gcd(f1c[1], f1c[2]) == Nemo.ZZ(1)

   @test divexact(f1c[2], f1c[1]) == Nemo.ZZ(3)

   @test f1c[2] - f1c[1] == Nemo.ZZ(2)

   @test f1c[1] + f1c[2] == Nemo.ZZ(4)

   @test length(string((x+y)^2)) > 3

   @test hash((x+y)^2) == hash(x^2+2*x*y+y^2)

   @test deepcopy(f1c[1]) == f1c[1]
end

@testset "Nemo.fq_nmod..." begin
   F, a = Nemo.FiniteField(7, 2, "a")

   R, (x, y) = PolynomialRing(F, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coeffs(f1)]
   @test isa(f1c[1], Singular.n_unknown{Nemo.fq_nmod})

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

   @test inv(f1c[2]) == F(5)

   @test gcd(f1c[1], f1c[2]) == F(1)

   @test divexact(f1c[2], f1c[1]) == F(3)

   @test f1c[2] - f1c[1] == F(2)

   @test f1c[1] + f1c[2] == F(4)

   @test length(string((x+a*y)^2)) > 3

   @test hash((x+a*y)^2) == hash(x^2+2*a*x*y+(a+4)*y^2)

   @test deepcopy(f1c[1]) == f1c[1]
end

@testset "Nemo.fq..." begin
   F, a = Nemo.FiniteField(Nemo.ZZ(7), 2, "a")

   R, (x, y) = PolynomialRing(F, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coeffs(f1)]
   @test isa(f1c[1], Singular.n_unknown{Nemo.fq})

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

   @test inv(f1c[2]) == F(5)

   @test gcd(f1c[1], f1c[2]) == F(1)

   @test divexact(f1c[2], f1c[1]) == F(3)

   @test f1c[2] - f1c[1] == F(2)

   @test f1c[1] + f1c[2] == F(4)

   @test length(string((x+a*y)^2)) > 3

   @test hash((x+a*y)^2) == hash(x^2+2*a*x*y+(a+4)*y^2)

   @test deepcopy(f1c[1]) == f1c[1]
end

@testset "Nemo.nf_elem..." begin
   U, z = Nemo.PolynomialRing(Nemo.QQ, "z")
   K, a = Nemo.NumberField(z^3 + 3z + 1, "a")

   R, (x, y) = PolynomialRing(K, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coeffs(f1)]
   @test isa(f1c[1], Singular.n_unknown{Nemo.nf_elem})

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

   @test inv(f1c[2]) == K(1)//3

   @test gcd(f1c[1], f1c[2]) == K(1)

   @test divexact(f1c[2], f1c[1]) == K(3)

   @test f1c[2] - f1c[1] == K(2)

   @test f1c[1] + f1c[2] == K(4)

   @test length(string((x+a*y)^2)) > 3

   @test hash((x+a*y)^2) == hash(x^2+2*a*x*y+a^2*y^2)

   @test deepcopy(f1c[1]) == f1c[1]
end

@testset "Nemo.NemoField..." begin
   U = Nemo.Generic.FractionField(Nemo.ZZ)

   R, (x, y) = PolynomialRing(U, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coeffs(f1)]
   @test isa(f1c[1], Singular.n_unknown{AbstractAlgebra.Generic.Frac{Nemo.fmpz}})

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

   @test inv(f1c[2]) == U(1, 3)

   @test gcd(f1c[1], f1c[2]) == U(1)

   @test divexact(f1c[2], f1c[1]) == U(3)

   @test f1c[2] - f1c[1] == U(2)

   @test f1c[1] + f1c[2] == U(4)

   @test length(string((x+y)^2)) > 3

   @test hash((x+y)^2) == hash(x^2+2*x*y+y^2)

   @test deepcopy(f1c[1]) == f1c[1]
end

@testset "Nemo.NemoRing..." begin
   U, z = Nemo.PolynomialRing(Nemo.ZZ, "z")

   R, (x, y) = PolynomialRing(U, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coeffs(f1)]
   @test isa(f1c[1], Singular.n_unknown{Nemo.fmpz_poly})

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

   @test gcd(f1c[1], f1c[2]) == U(1)

   @test divexact(f1c[2], f1c[1]) == U(3)

   @test f1c[2] - f1c[1] == U(2)

   @test f1c[1] + f1c[2] == U(4)

   @test length(string((x+z*y)^2)) > 3

   @test hash((x+z*y)^2) == hash(x^2+2*z*x*y+z^2*y^2)

   @test deepcopy(f1c[1]) == f1c[1]
end
