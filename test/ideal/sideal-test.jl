function test_sideal_constructors()
   print("sideal.constructors...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x, y)
   S = parent(I)

   @test elem_type(S) == sideal{spoly{n_Q}}
   @test elem_type(IdealSet{spoly{n_Q}}) == sideal{spoly{n_Q}}
   @test parent_type(sideal{spoly{n_Q}}) == IdealSet{spoly{n_Q}}

   @test base_ring(S) == R
   @test base_ring(I) == R

   @test typeof(S) <: AbstractAlgebra.Set
   @test typeof(I) <: AbstractAlgebra.Module

   I1 = Ideal(R)
   I2 = Ideal(R, x + y)
   I3 = Ideal(R, x + y, y^2 + 2)
   I4 = Ideal(R, [x + y, y^2 + 2])
   I5 = MaximalIdeal(R, 0)
   I6 = MaximalIdeal(R, 4)

   @test isa(I1, sideal)
   @test isa(I2, sideal)
   @test isa(I3, sideal)
   @test isa(I4, sideal)
   @test isa(I5, sideal)
   @test isa(I6, sideal)

   println("PASS")
end

function test_sideal_manipulation()
   print("sideal.manipulation...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I0 = Ideal(R)
   I1 = Ideal(R, x, y)

   @test ngens(I0) == 0
   @test ngens(I1) == 2
   @test x in gens(I1) && y in gens(I1)

   I2 = deepcopy(I1)

   @test isequal(I1, I2)

   @test I2[2] == y
   I2[2] = x^2 + 1
   @test I2[2] == x^2 + 1

   @test iszero(I0)

   @test iszerodim(I1)

   @test isconstant(Ideal(R, R(1), R(2)))

   @test isvar_generated(Ideal(R, x))
   @test isvar_generated(Ideal(R, y))
   @test isvar_generated(Ideal(R, x, y))
   @test !isvar_generated(Ideal(R, R(1)))
   @test !isvar_generated(Ideal(R, x + y))

   println("PASS")
end

function test_sideal_binary_ops()
   print("sideal.binary_ops...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I1 = Ideal(R, x)
   I2 = Ideal(R, y)
   I3 = Ideal(R, x, y)
   I4 = Ideal(R, x*y)

   @test equal(I1 + I2, I3)
   @test equal(I1*I2, I4)

   println("PASS")
end

function test_sideal_powering()
   print("sideal.powering...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x^2, x*y + 1)

   @test equal(I^0, Ideal(R, R(1)))

   S = I

   for i = 1:5
      @test equal(S, I^i)
      S *= I
   end

   println("PASS")
end

function test_sideal_containment()
   print("sideal.containment...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   @test contains(Ideal(R, x, y), Ideal(R, x))
   @test !contains(Ideal(R, x), Ideal(R, x, y))

   println("PASS")
end

function test_sideal_comparison()
   print("sideal.comparison...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   @test equal(Ideal(R, x, y, x - y), Ideal(R, x, y, x + y))
   @test !equal(Ideal(R, x), Ideal(R, x, y))
   @test !isequal(Ideal(R, x, y, x - y), Ideal(R, x, y, x + y))
   @test isequal(Ideal(R, x, y), Ideal(R, x, y))

   println("PASS")
end

function test_sideal_leading_terms()
   print("sideal.leading_terms...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   @test equal(lead(Ideal(R, x^2 + x*y + 1, 2y^2 + 3, R(7))), Ideal(R, x^2, 2y^2, R(7)))

   println("PASS")
end

function test_sideal_local()
   print("sideal.local...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:negdegrevlex)

   I = Ideal(R, y, x^2, (1 + y^3) * (x^2 - y))
   M = Singular.minimal_generating_set(I)

   @test equal(I, Ideal(R, x^2, y))
   println("PASS")
end

function test_sideal_intersection()
   print("sideal.intersection...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I1 = Ideal(R, x^2 + x*y + 1, 2y^2 + 3, R(7))
   I2 = Ideal(R, x*y^2 + x + 1, 2x*y + 1, 7x + 1)

   I = intersection(I1, I2)

   @test contains(I1, I)
   @test contains(I2, I)

   println("PASS")
end

function test_sideal_quotient()
   print("sideal.quotient...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x^2 + x*y + 1, 2y^2 + 3)
   J = Ideal(R, x*y^2 + x + 1, 2x*y + 1)
   K = Ideal(R, x*y + 1)

   A = quotient(I, J + K)
   B = intersection(quotient(I, J), quotient(I, K))

   @test equal(A, B)

   println("PASS")
end

 function test_sideal_saturation()
   print("sideal.saturation...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, (x^2 + x*y + 1)*(2y^2+1)^3, (2y^2 + 3)*(2y^2+1)^2)
   J = Ideal(R, 2y^2 + 1)

   @test equal(saturation(I, J), Ideal(R, 2y^2 + 3, x^2 + x*y + 1))

   I = Ideal(R, (x*y + 1)*(2x^2*y^2 + x*y - 2) + 2x*y^2 + x, 2x*y + 1)
   J = Ideal(R, x)

   @test equal(satstd(I, J), std(saturation(I, J)))

   println("PASS")
end

function test_sideal_slimgb()
   print("sideal.slimgb...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x^2 + x*y + 1, 2y^2 + 3)
   J = Ideal(R, 2*y^2 + 3, x^2 + x*y + 1)

   A = slimgb(I)

   @test isequal(lead(A), Ideal(R, 2y^2, x^2)) ||
         isequal(lead(A), Ideal(R, x^2, 2*y^2))
   @test A.isGB == true

   B = slimgb(I, complete_reduction=true)

   @test isequal(B, Ideal(R, 2y^2 + 3, x^2 + x*y + 1)) ||
         isequal(B, Ideal(x^2 + x*y + 1, 2y^2 + 3))
   @test B.isGB == true

   println("PASS")
end

function test_sideal_std()
   print("sideal.std...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x^2 + x*y + 1, 2y^2 + 3)
   J = Ideal(R, 2*y^2 + 3, x^2 + x*y + 1)

   A = std(I)

   @test isequal(lead(A), Ideal(R, 2y^2, x^2)) ||
         isequal(lead(A), Ideal(R, x^2, 2*y^2))
   @test A.isGB == true

   B = std(I, complete_reduction=true)

   @test isequal(B, Ideal(R, 2y^2 + 3, x^2 + x*y + 1)) ||
         isequal(B, Ideal(x^2 + x*y + 1, 2y^2 + 3))
   @test B.isGB == true

   println("PASS")
end

function test_sideal_reduction()
   print("sideal.reduction...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   f = x^2*y + 2y + 1
   g = y^2 + 1

   I = Ideal(R, (x^2 + 1)*f + (x + y)*g + x + 1, (2y^2 + x)*f + y)
   J = std(Ideal(R, f, g))

   @test isequal(reduce(I, J), Ideal(R, x + 1, y))

   h = (x^2 + 1)*f + (x + y)*g + x + 1

   @test reduce(h, J) == x + 1

   println("PASS")
end

function test_sideal_free_resolution()
   print("sideal.free_resolution...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x^2*y + 2y + 1, y^2 + 1)

   F1 = fres(std(I), 4)
   F2 = sres(std(I), 4)

   # check resolution is of the correct length
   @test length(F1) == 2
   @test length(F2) == 2

   M1 = Singular.Matrix(F1[1])
   N1 = Singular.Matrix(F1[2])

   M2 = Singular.Matrix(F2[1])
   N2 = Singular.Matrix(F2[2])

   # check we have a complex
   @test iszero(M1*N1)
   @test iszero(M2*N2)

   println("PASS")
end

function test_sideal_syzygy()
   print("sideal.syzygy...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x^2*y + 2y + 1, y^2 + 1)

   F = syz(I)

   M = Singular.Matrix(I)
   N = Singular.Matrix(F)

   # check they are actually syzygies
   @test iszero(M*N)

   println("PASS")
end

function test_sideal_kernel()
   print("sideal.kernel...")

   # twisted cubic

   P1, (t_0, t_1) = PolynomialRing(QQ, ["t_0", "t_1"])
   P3, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"])
   I = Ideal(P1, t_0^3, t_0^2*t_1, t_0*t_1^2, t_1^3)
   J = kernel(P3, I)
   @test ngens(J) == 3 && J[1] == z^2-y*w && J[2] == y*z-x*w && J[3] == y^2-x*z
   println("PASS")
end

function test_sideal_eliminate()
   print("sideal.eliminate...")

   R, (x, y, t) = PolynomialRing(QQ, ["x", "y", "t"])

   I = Ideal(R, x - t^2, y - t^3)

   J = eliminate(I, t)

   @test equal(J, Ideal(R, x^3 - y^2))

   K = eliminate(I, t, x)

   @test equal(K, Ideal(R))

   println("PASS")
end

function test_sideal_jet()
   print("sideal.jet...")

   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

   I = Ideal(R, x^5 - y^2, y^3 - x^6 + z^3)

   J = jet(I, 3)

   @test equal(J, Ideal(R, - y^2, y^3 + z^3))

   println("PASS")
end

function test_sideal_jacobi()
   print("sideal.jacob...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, 2*x + 5*y, 2*y + 5*x)

   J = jacobi(I)

   M = R[2 5; 5 2]

   @test J == M

   println("PASS")
end

function test_sideal_zerodim()
   print("sideal.zerodim...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"]; ordering=:negdegrevlex)

   I = Ideal(R, 3*x^2 + y^3, x*y^2)

   I = std(I)

   dim = vdim(I)

   J = kbase(I)

   B = Ideal(R, R(1), x, y, x*y, y^2, y^3, y^4)

   f = highcorner(I)

   # Check dimension
   @test dim == 7

   # Check vector space basis
   @test equal(J, B)

   #Check highcorner
   @test f == y^4

   println("PASS")
end

function test_sideal_extend()
   print("sideal.extend...")

   R1 , (x, y, z) = Singular.PolynomialRing(Singular.QQ,
   ["x", "y", "z"]; ordering=:negdegrevlex)

   I1 = Ideal(R1, x^2+y, z*y+5)

   S1, (x, y, z) = Singular.PolynomialRing(Nemo.QQ,
   ["x", "y", "z"]; ordering=:negdegrevlex)

   I2 = Ideal(S1, x^2+y, z*y+5)

   I3 = extend(I1, Nemo.QQ, S1)

   R2 , (x, y, z) = Singular.PolynomialRing(Singular.ZZ,
   ["x", "y", "z"]; ordering=:negdegrevlex)

   I4 = Ideal(R2, x^2+y, z*y+5)

   S2, (x, y, z) = Singular.PolynomialRing(Singular.Fp(5),
   ["x", "y", "z"]; ordering=:negdegrevlex)

   I5 = Ideal(S2, x^2+y, z*y+5)

   I6 = extend(I4, Singular.Fp(5), S2)

   @test equal(I2, I3)

   @test equal(I5, I6)


   println("PASS")
end

function test_sideal()
   test_sideal_constructors()
   test_sideal_manipulation()
   test_sideal_binary_ops()
   test_sideal_powering()
   test_sideal_containment()
   test_sideal_comparison()
   test_sideal_local()
   test_sideal_leading_terms()
   test_sideal_intersection()
   test_sideal_quotient()
   test_sideal_saturation()
   test_sideal_slimgb()
   test_sideal_std()
   test_sideal_reduction()
   test_sideal_free_resolution()
   test_sideal_syzygy()
   test_sideal_eliminate()
   test_sideal_kernel()
   test_sideal_jet()
   test_sideal_jacobi()
   test_sideal_zerodim()
   test_sideal_extend()

   println("")
end

