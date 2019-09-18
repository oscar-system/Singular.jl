@testset "spoly.constructors..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   @test elem_type(R) == spoly{n_Z}
   @test elem_type(PolyRing{n_Z}) == spoly{n_Z}
   @test parent_type(spoly{n_Z}) == PolyRing{n_Z}
   @test base_ring(R) == ZZ

   typeof(R) <: Nemo.Ring

   a = R()

   @test base_ring(a) == ZZ
   @test parent(a) == R

   @test isa(a, spoly)

   b = R(123)

   @test isa(b, spoly)

   c = R(BigInt(123))

   @test isa(c, spoly)

   d = R(c)

   @test isa(d, spoly)

   f = R(Nemo.ZZ(123))

   @test isa(f, spoly)

   g = R(ZZ(123))

   @test isa(g, spoly)

   S, (y, ) = PolynomialRing(QQ, ["y", ])

   h = S(ZZ(123))

   @test isa(h, spoly)

   T, (z, ) = PolynomialRing(Nemo.ZZ, ["z", ])

   k = T(123)

   @test isa(k, spoly)

   S = @PolynomialRing(ZZ, "x", 50)

   @test isa(x17, spoly)

   T = @PolynomialRing(ZZ, "y", 50, :lex)

   @test isa(y7, spoly)

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   @test isgen(x)
   @test isgen(y)
   @test !isgen(R(1))
   @test !isgen(R(0))
   @test !isgen(2x)
   @test !isgen(x + y)
   @test !isgen(x^2)

   @test has_global_ordering(R)
   @test !has_local_ordering(R)
   @test !has_mixed_ordering(R)

   @test length(symbols(R)) == 2
   @test symbols(R) == [:x, :y]

   R, (x, y) = PolynomialRing(ZZ, ["x", "y"]; ordering=:lex)

   M = MPolyBuildCtx(R)
   push_term!(M, ZZ(2), [1, 2])
   push_term!(M, ZZ(1), [1, 1])
   push_term!(M, ZZ(2), [3, 2])
   f = finish(M)

   @test f == 2*x^3*y^2+2*x*y^2+x*y
end

@testset "spoly.printing..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   @test string(3x^2 + 2x + 1) == "3*x^2+2*x+1"
end

@testset "spoly.manipulation..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   @test isone(one(R))
   @test iszero(zero(R))
   @test isunit(R(1)) && isunit(R(-1))
   @test !isunit(R(2)) && !isunit(R(0)) && !isunit(x)
   @test isgen(x)
   @test !isgen(R(1)) && !isgen(x + 1)
   @test isconstant(R(0)) && isconstant(R(1))
   @test !isconstant(x) && !isconstant(x + 1)
   @test ismonomial(x) && ismonomial(R(1))
   @test !ismonomial(2x) && !ismonomial(x + 1)
   @test isterm(2x) && !isterm(x + 1)
   @test length(x^2 + 2x + 1) == 3
   @test total_degree(x^2 + 2x + 1) == 2
   @test order(x^2 + 2x + 1) == 0

   @test lead_exponent(x^3 + 2x + 1) == [3]

   @test deepcopy(x + 2) == x + 2

   @test characteristic(R) == 0

   @test nvars(R) == 1
   pol = x^5 + 3x + 2

   @test length(collect(coeffs(pol))) == length(pol)
   @test length(collect(exponent_vectors(pol))) == length(pol)

   polzip = zip(coeffs(pol), monomials(pol), terms(pol))
   r = R()
   for (c, m, t) in polzip
      r += c*m
      @test t == c*m
   end

   @test pol == r

   R, (x, ) = PolynomialRing(ResidueRing(ZZ, 6), ["x", ])

   @test characteristic(R) == 6

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])
   p = x + y
   q = x
   @test Singular.substitute_variable(p, 2, q) == 2*x
   @test Singular.permute_variables(q, [2, 1], R) == y

   for i = 1:nvars(R)
      @test gen(R, i) == gens(R)[i]
   end
   @test gen(R, 1) == x

   @test ordering(R) == :degrevlex
   @test degree(x^2*y^3 + 1, 1) == 2
   @test degree(x^2*y^3 + 1, y) == 3
   @test degree(R(), 1) == -1
   @test degrees(x^2*y^3) == [2, 3]
   @test vars(x^2 + 3x + 1) == [x]
   @test var_index(x) == 1 && var_index(y) == 2
   @test lc(3x^2 + 2x + 1) == 3
   @test lm(3x^2 + 2x + 1) == x^2
   @test lt(3x^2 + 2x + 1) == 3x^2
end

@testset "spoly.change_base_ring..." begin
   R1, (x, ) = PolynomialRing(ZZ, ["x", ])

   a1 = x^2 + 3x + 1

   b1 = change_base_ring(a1, QQ)

   @test isa(b1, spoly{n_Q})

   R2, (x, y) = PolynomialRing(ZZ, ["x", "y"])
   a2 = x^5+y^3+1

   R3, (x, y) = PolynomialRing(QQ, ["x", "y"])
   a3 = x^5+y^3+1

   R4, (x, y) = PolynomialRing(Nemo.QQ, ["x", "y"])
   a4 = x^5+y^3+1

   a5 = change_base_ring(a2, QQ)
   a6 = change_base_ring(a2, CoefficientRing(Nemo.QQ))

   @test a3 == a5
   @test a4 == a6
end

@testset "spoly.multivariate_coeff..." begin
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

   f = 2x^2*y^2 + 3x*y^2 - x^2*y + 4x*y - 5y + 1

   @test coeff(f, [2], [1]) == -x^2 + 4x - 5
   @test coeff(f, [y], [1]) == -x^2 + 4x - 5
end

@testset "spoly.unary_ops..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1

   @test -a == -x^2 - 3x - 1
end

@testset "spoly.binary_ops..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1
   b = 2x + 4

   @test a + b == x^2+5*x+5
   @test a - b == x^2+x-3
   @test a*b == 2*x^3+10*x^2+14*x+4
end

@testset "spoly.comparison..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1

   @test a == deepcopy(a)
   @test a != x
end

@testset "spoly.powering..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1

   @test a^0 == 1
   @test a^1 == x^2 + 3x + 1
   @test a^3 == x^6+9*x^5+30*x^4+45*x^3+30*x^2+9*x+1
end

@testset "spoly.exact_division..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1
   b = 2x + 4

   @test divexact(a*b, a) == b
end

@testset "spoly.adhoc_exact_division..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   a = x^2 + 3x + 1

   @test divexact(2a, 2) == a
   @test divexact(2a, BigInt(2)) == a
   @test divexact(2a, ZZ(2)) == a

   R, (x, ) = PolynomialRing(QQ, ["x", ])

   a = x^2 + 3x + 1

   @test divexact(2a, 2) == a
   @test divexact(2a, BigInt(2)) == a
   @test divexact(2a, ZZ(2)) == a
   @test divexact(2a, 2//3) == 3a
   @test divexact(2a, QQ(2//3)) == 3a
   @test divexact(2a, BigInt(2)//3) == 3a
end

@testset "spoly.euclidean_division..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   a = x^2*y^2 + 3x + 1
   b = x*y + 1

   q, r = divrem(a, b)
   @test a == b*q + r

   q2 = div(a, b)
   @test q2 == q
end

@testset "spoly.divides..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   a = x^2 + 3x + 1
   b = x*y + 1

   flag, q = divides(a*b, b)
   @test flag && q == a

   flag, q = divides(a, y)
   @test !flag

   val, q = remove(a*b^3, b)
   @test val == 3 && q == a
   @test valuation(a*b^3, b) == 3
end

@testset "spoly.gcd_lcm..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])
   a = x^2 + 3x + 1
   b = 2x + 4
   c = 2x^2 + 1

   @test gcd(a*c, b*c) == c

   @test lcm(a, b) == a*b

   @test primpart(2*a) == a

   @test content(2*a) == 2
end

@testset "spoly.extended_gcd..." begin
   R, (x, ) = PolynomialRing(QQ, ["x", ])

   a = x^2 + 3x + 1
   b = 2x + 4

   g, s, t = gcdx(a, b)

   @test s*a + t*b == g
end

@testset "spoly.evaluate..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   f = x^2*y + 2x + 1

   @test evaluate(f, [2, 3]) == 17
   @test evaluate(f, [ZZ(2), ZZ(3)]) == 17
   @test evaluate(f, [QQ(2), QQ(3)]) == 17
   @test evaluate(f, [x + 1, y - 1]) == x^2*y - x^2 + 2*x*y + y + 2
   @test evaluate(f, [x], [1]) == y + 3
   @test f(2, 3) == 17
end

@testset "spoly.inflation_deflation..." begin
   R, (x, y) = PolynomialRing(ZZ, ["x", "y"])

   f = x^7*y^7 + 3x^4*y^4 + 2x*y

   @test inflate(deflate(f, deflation(f)...), deflation(f)...) == f
end

@testset "spoly.Polynomials..." begin
   R, (x, ) = PolynomialRing(ZZ, ["x", ])

   S, y = Nemo.PolynomialRing(R, "y")

   f = (1 + 2x + 3x^2)*y + (2x + 3)

   g = f^2

   @test g == (9*x^4+12*x^3+10*x^2+4*x+1)*y^2+(12*x^3+26*x^2+16*x+6)*y+(4*x^2+12*x+9)
end

@testset "spoly.convert_MPoly_to_SingularPoly..." begin
   for num_vars = 2:10
      var_names = ["x$j" for j in 1:num_vars]
      ord = AbstractAlgebra.rand_ordering()

      R, vars_R = AbstractAlgebra.Generic.PolynomialRing(Nemo.ZZ, var_names; ordering=ord)
      Rsing, vars_Rsing = Singular.AsEquivalentSingularPolynomialRing(R)
      for iter in 1:10
         f = AbstractAlgebra.Generic.rand(R, 5:10, 1:10, -100:100)
         g = AbstractAlgebra.Generic.rand(R, 5:10, 1:10, -100:100)
         @test Rsing(f * g) == Rsing(f) * Rsing(g)
         @test Rsing(f + g) == Rsing(f) + Rsing(g)
         @test Rsing(f - g) == Rsing(f) - Rsing(g)
         @test R(Rsing(f)) == f
      end

      S, vars_S = Singular.PolynomialRing(Nemo.ZZ, var_names; ordering=ord)
      SAA, vars_SAA = Singular.AsEquivalentAbstractAlgebraPolynomialRing(S)
      # FIXME: Better generate random polynomials
      f = 1 + vars_S[1]*vars_S[2]
      g = vars_S[1]^2 + vars_S[2]^3
      @test SAA(f * g) == SAA(f) * SAA(g)
      @test SAA(f + g) == SAA(f) + SAA(g)
      @test SAA(f - g) == SAA(f) - SAA(g)
      @test S(SAA(f)) == f
   end
end

@testset "spoly.test_spoly_differential..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   f = x^3 + y^6

   J = jacobi(f)
   f1 = derivative(f, 1)
   f2 = derivative(f, y)
   jf = jet(f, 3)

   I = Ideal(R, x^2, y^5)

   # Check derivative
   @test f1 == 3*x^2
   @test f2 == 6*y^5

   #Check jacobi
   @test equal(I, J)

   #Check jet
   @test jf == x^3
end

@testset "spoly.test_spoly_factor..." begin
   R1 , (x, y, z, w) = Singular.PolynomialRing(Singular.QQ,
   ["x", "y", "z", "w"]; ordering=:negdegrevlex)
   f1 = 113*(2*y^7 + w^2)^3*(1 + x)^2*(x + y*z)^2

   R2 , (x, y, z, w) = Singular.PolynomialRing(Singular.ZZ,
   ["x", "y", "z", "w"]; ordering=:negdegrevlex)
   f2 = 123*(57*y^3 + w^5)^3*(x^2 + x+1)^2*(x + y*z)^2

   R3 , (x, y, z, w) = Singular.PolynomialRing(Singular.Fp(3),
   ["x", "y", "z", "w"]; ordering=:negdegrevlex)
   f3 = 7*(y^3 + w^3)*(1 + x)^2*(x + y*z)^2

   for f in [f1, f2, f3]
      #Test factor
      F = factor(f)
      g = F.unit
      for p in keys(F.fac)
         g = g*p^F[p]
      end
     @test f == g

     #Test factor_squarefree over QQ and Fp
      if typeof(f.parent) != PolyRing{n_Z}
         F = factor_squarefree(f)
         g = F.unit
         for p in keys(F.fac)
            g = g*p^F[p]
         end
         @test f == g
      end
   end
end

@testset "spoly.hash..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   @test hash(x) == hash(x+y-y)
   @test hash(x,zero(UInt)) == hash(x+y-y,zero(UInt))
   @test hash(x,one(UInt)) == hash(x+y-y,one(UInt))
end

