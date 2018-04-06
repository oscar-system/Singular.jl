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
   test_sideal_kernel()

   println("")
end

