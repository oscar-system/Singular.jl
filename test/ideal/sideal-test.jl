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

   I2 = deepcopy(I1)

   @test isequal(I1, I2)

   @test I2[2] == y
   I2[2] = x^2 + 1
   @test I2[2] == x^2 + 1

   @test iszero(I0)

   @test iszerodim(I1)

   @test isconstant(Ideal(R, R(1), R(2)))

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

function test_sideal()
   test_sideal_constructors()
   test_sideal_manipulation()
   test_sideal_binary_ops()
   test_sideal_kernel()

   println("")
end

