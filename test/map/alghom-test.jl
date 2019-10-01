function test_alghom_constructors()
   print("alghom.constructors...")
   L = Singular.FiniteField(3, 2, String("a"))
   R, (x, y, z, w) = Singular.PolynomialRing(L[1], ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)
   S, (a, b, c) = Singular.PolynomialRing(L[1], ["a", "b", "c"];
                             ordering=:degrevlex)
   I = Ideal(S, [a, a + b^2, b - c, c + b])
   f = Singular.AlgebraHomomorphism(R, S, gens(I), I.ptr)

   @test f.domain == R
   @test f.codomain == S
   @test isequal(f.image, gens(I))

   f = Singular.AlgebraHomomorphism(R, S, gens(I))
   @test f.domain == R
   @test f.codomain == S
   @test isequal(f.image, gens(I))

   println("PASS")
end

function test_alghom_apply()
   print("alghom.apply...")
   R, (x, y, z, w) = Singular.PolynomialRing(Singular.Fp(3), ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)
   S, (a, b, c) = Singular.PolynomialRing(Singular.Fp(3), ["a", "b", "c"];
                             ordering=:degrevlex)
   I = Ideal(S, [a, a + b^2, b - c, c + b])
   f = Singular.AlgebraHomomorphism(R, S, gens(I), I.ptr)
   id  = Singular.IdentityAlgebraHomomorphism(S)

   J = Ideal(R, [x, y^3])
   p = x + y^3 + z*w

   J1 = Ideal(S, [a, a^3 + b^6])
   p1 = b^6 + a^3 + b^2 - c^2 + a

   @test isequal(f(J), J1)
   @test f(p) == p1
   @test isequal(gens(id(I)), gens(I))
   @test id(f(p)) == f(p)

   println("PASS")
end

function test_alghom_compose()
   print("alghom.compose...")
   R, (x, y, z, w) = Singular.PolynomialRing(Singular.QQ, ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)
   S, (a, b, c) = Singular.PolynomialRing(Singular.QQ, ["a", "b", "c"];
                             ordering=:degrevlex)
   V = [a, a + b^2, b - c, c + b]
   W = [x^2, x + y + z, z*y]

   L = [x^2, 2*x^2 + 2*x*y + y^2 + 2*x*z + 2*y*z + z^2, x + y + z - y*z,
                x + y + z + y*z]

   f = Singular.AlgebraHomomorphism(R, S, V)
   g = Singular.AlgebraHomomorphism(S, R, W)
   idR  = Singular.IdentityAlgebraHomomorphism(R)
   idS  = Singular.IdentityAlgebraHomomorphism(S)

   h1 = f*g
   h2 = f*idS
   h3 = idR*f
   h4 = idR*idR

   @test h1.domain == R && h1.codomain == R && isequal(h1.image, L)
   @test h2.domain == R && h2.codomain == S && isequal(h2.image, V)
   @test h3.domain == R && h3.codomain == S && isequal(h3.image, V)
   @test h4.domain == R && isequal(h4.image, gens(MaximalIdeal(R, 1)))

   println("PASS")
end

function test_alghom_preimage()
   print("alghom.preimage...")
   R, (x, y, z, w) = Singular.PolynomialRing(Singular.QQ, ["x", "y", "z", "w"];
                             ordering=:negdegrevlex)
   S, (a, b, c) = Singular.PolynomialRing(Singular.QQ, ["a", "b", "c"];
                             ordering=:degrevlex)
   I = Ideal(S, [a, a + b^2, b - c, c + b])

   f = Singular.AlgebraHomomorphism(R, S, gens(I), I.ptr)
   idS  = Singular.IdentityAlgebraHomomorphism(S)

   @test isequal(Singular.preimage(f, I), MaximalIdeal(R, 1))
   @test isequal(Singular.kernel(f), Ideal(R, 4*x - 4*y + z^2 + 2*z*w + w^2))
   @test isequal(Singular.preimage(idS, I), I)
   @test isequal(Singular.kernel(idS), Ideal(S,))

   println("PASS")
end

function test_alghom()
   test_alghom_constructors()
   test_alghom_apply()
   test_alghom_compose()
   test_alghom_preimage()

   println("")
end
