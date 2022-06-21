@testset "quotient.PolynomialRing" begin
   R, (x, y, z) = @inferred PolynomialRing(QQ, ["x", "y", "z"])
   Q, (x, y, z) = @inferred QuotientRing(R, Ideal(R, x^2 - y*z))

   @test is_quotient_ring(Q)

   I = Ideal(Q, x)
   @test !I.isGB
   @test iszero(reduce(Ideal(Q, x, y*z), std(I)))
end

@testset "quotient.GAlgebra" begin

   R, _ = @inferred PolynomialRing(QQ, ["z", "u", "v", "w"])

   G, (z, u, v, w) = @inferred GAlgebra(R, Singular.Matrix(R, [0 -1 -1 -1; 0 0 -1 -1; 0 0 0 -1; 0 0 0 0]),
                                           Singular.Matrix(R, [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]))

   # error thrown if ideal is not a Groebner basis (and not two sided)
   # note principal ideals are now auto recognized as isGB=true
   I = Ideal(G, [u^2-v^3*w, z^2*v-u^2])
   @test_throws ErrorException QuotientRing(G, I)

   # error thrown if ideal is GB but not two sided
   I = std(I)
   @test_throws ErrorException QuotientRing(G, I)

   # error thrown if ideal is two sided but not a GB
   I = Ideal(G, [u^2-v^3*w, z^2*v-u^2], twosided = true)
   @test_throws ErrorException QuotientRing(G, I)

   # all good if two sided and GB
   I = std(I)
   @test QuotientRing(G, I)[1] isa PluralRing

   # another test
   I = @inferred Ideal(G, [z^2, u^2, v^2, w^2, z*u*v - w]; twosided = true)
   Q, (z, u, v, w) = @inferred QuotientRing(G, std(I))
   @test Q isa PluralRing
end

@testset "quotient.FreeAlgebra" begin
   R, (x, d) = FreeAlgebra(Fp(3), ["x", "d"], 5)
   I = @inferred Ideal(R, [x^4, d^3, d*x - x*d - 1])
   @test length(gens(std(I))) == 3

   I = @inferred Ideal(R, [d*x - x*d - 1])
   Q, (x, d) = @inferred QuotientRing(R, std(I))
   @test length(@inferred gens(std(Ideal(Q, [x^4, d^3])))) == 2
   @test length(@inferred gens(std(Ideal(Q, [x^2, d^3])))) == 1
end

