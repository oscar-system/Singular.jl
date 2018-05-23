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

function test_nemo()
   test_nemo_fmpq()

   println("")
end
