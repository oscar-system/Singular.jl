@testset "Singular.LibPuiseuxexpansions" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   f = y^3 + x^2 + x^8

   L = Singular.LibPuiseuxexpansions.puiseux(f,8,1)

   @test L[1][1] == 1//9*x^38 - 1//3*x^20 - x^2
   @test L[1][2] == 3
   @test L[1][3] == nothing
end
