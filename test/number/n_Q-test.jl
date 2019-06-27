function test_n_Q_constructors()
   print("n_Q.constructors...")

   @test elem_type(QQ) == n_Q
   @test elem_type(Rationals) == n_Q
   @test parent_type(n_Q) == Rationals
   @test base_ring(QQ) == ZZ

   typeof(QQ) <: Nemo.Field

   a = QQ()

   @test base_ring(a) == ZZ
   @test parent(a) == QQ

   @test isa(a, n_Q)

   b = QQ(123)

   @test isa(b, n_Q)

   c = QQ(BigInt(123))

   @test isa(c, n_Q)

   d = QQ(c)

   @test isa(d, n_Q)

   f = QQ(Nemo.ZZ(123))

   @test isa(f, n_Q)

   g = QQ(ZZ(123))

   @test isa(g, n_Q)

   h = QQ(1, 2)

   @test isa(h, n_Q)

   i = QQ(1//2)

   @test isa(i, n_Q)
   @test i == QQ(1) // QQ(2)
   
   j = QQ(BigInt(1)//BigInt(2))

   @test isa(j, n_Q)
   @test j == QQ(1) // QQ(2)

   println("PASS")
end

function test_n_Q_printing()
   print("n_Q.printing...")

   @test string(QQ(123)) == "123"

   println("PASS")
end

function test_n_Q_manipulation()
   print("n_Q.manipulation...")

   @test isone(one(QQ))
   @test iszero(zero(QQ))
   @test isunit(QQ(1)) && isunit(QQ(2))
   @test !isunit(QQ(0))

   @test numerator(QQ(2)) == 2
   @test denominator(QQ(2)) == 1
   @test denominator(QQ()) == 1

   @test abs(QQ(-2)) == 2
   @test abs(QQ(2)) == 2
   @test abs(QQ()) == 0

   @test deepcopy(QQ(2)) == QQ(2)

   println("PASS")
end

function test_n_Q_unary_ops()
   print("n_Q.unary_ops...")

   @test -QQ(123) == -123
   @test -QQ() == 0

   println("PASS")
end

function test_n_Q_binary_ops()
   print("n_Q.binary_ops...")

   a = QQ(2)
   b = QQ(3)

   @test a + b == 5
   @test a - b == -1
   @test a*b == 6

   println("PASS")
end

function test_n_Q_adhoc_binary()
   print("n_Q.adhoc_binary...")

   @test QQ(2) + 3 == 5
   @test 2 + QQ(3) == 5
   @test QQ(2) - 3 == -1
   @test 2 - QQ(3) == -1
   @test 2*QQ(3) == 6
   @test QQ(2)*3 == 6

   println("PASS")
end

function test_n_Q_comparison()
   print("n_Q.comparison...")

   @test QQ(2) < QQ(3)
   @test QQ(2) <= QQ(3)
   @test QQ(3) > QQ(2)
   @test QQ(3) >= QQ(3)
   @test QQ(2) == QQ(2)
   @test isequal(QQ(2), QQ(2))

   println("PASS")
end

function test_n_Q_adhoc_comparison()
   print("n_Q.ad_hoc_comparison...")

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

   println("PASS")
end

function test_n_Q_powering()
   print("n_Q.powering...")

   @test QQ(2)^10 == 1024

   println("PASS")
end

function test_n_Q_exact_division()
   print("n_Q.exact_division...")

   @test divexact(QQ(12), QQ(2)) == 6

   println("PASS")
end

function test_n_Q_gcd_lcm()
   print("n_Q.gcd_lcm...")

   @test gcd(QQ(6), QQ(12)) == 1
   @test gcd(QQ(-6), QQ(12)) == 1
   @test gcd(QQ(6), QQ(-12)) == 1
   @test gcd(QQ(-6), QQ(-12)) == 1

   @test gcd(QQ(0), QQ(0)) == 0

   @test lcm(QQ(4), QQ(6)) == 24

   println("PASS")
end

function test_n_Q_reconstruct()
   print("n_Q.reconstruct...")

   @test reconstruct(ZZ(11), ZZ(15)) == divexact(QQ(-1), QQ(4))
   @test reconstruct(ZZ(11), 15) == divexact(QQ(-1), QQ(4))
   @test reconstruct(11, ZZ(15)) == divexact(QQ(-1), QQ(4))

   println("PASS")
end

function test_n_Q_Polynomials()
   print("n_Q.Polynomials...")

   R, x = Nemo.PolynomialRing(QQ, "x")

   f = 1 + 2x + 3x^2

   g = f^2

   @test g == 9*x^4+12*x^3+10*x^2+4*x+1

   println("PASS")
end

function test_n_Q_conversions()
   print("n_Q.conversions...")

   @test Singular.QQ(1//2) == Singular.QQ(1) // Singular.QQ(2)
   @test AbstractAlgebra.QQ(Singular.QQ(1) // Singular.QQ(2)) == 1//2

   println("PASS")
end

function test_n_Q()
   test_n_Q_constructors()
   test_n_Q_printing()
   test_n_Q_manipulation()
   test_n_Q_unary_ops()
   test_n_Q_binary_ops()
   test_n_Q_adhoc_binary()
   test_n_Q_comparison()
   test_n_Q_adhoc_comparison()
   test_n_Q_powering()
   test_n_Q_exact_division()
   test_n_Q_gcd_lcm()
   test_n_Q_reconstruct()
   test_n_Q_Polynomials()
   test_n_Q_conversions()

   println("")
end

