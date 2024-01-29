@testset "Nemo.QQFieldElem" begin
   R, (x, y) = polynomial_ring(Nemo.QQ, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]

# disable this until we figure out if QQ is mutable or not
#   @test isa(f1c[1], Singular.n_FieldElem{Singular.FieldElemWrapper{Nemo.QQField, Nemo.QQFieldElem}})
   @test isa(Nemo.QQ(f1c[1]), Nemo.QQFieldElem)
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]

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

@testset "Nemo.ZZRingElem" begin
   R, (x, y) = polynomial_ring(Nemo.ZZ, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
# disable this until we figure out if ZZ is mutable or not
#   @test f1c[1] isa Singular.n_RingElem{Singular.RingElemWrapper{Nemo.ZZRing, Nemo.ZZRingElem}}
   @test Nemo.ZZ(f1c[1]) isa Nemo.ZZRingElem
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]

   @test string(first(coefficients(f3))) == "1"

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

@testset "Nemo.fqPolyRepFieldElem" begin
   F, _ = Nemo.Native.finite_field(Nemo.next_prime(fld(typemax(Int),2)), 2, "a")
   R, _ = polynomial_ring(F, ["x", "y"])
   @test R isa Singular.PolyRing{Singular.n_FieldElem{Nemo.fqPolyRepFieldElem}}

   F, a = Nemo.Native.finite_field(7, 2, "a")

   R, (x, y) = polynomial_ring(F, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa Singular.n_FieldElem{Nemo.fqPolyRepFieldElem}
   @test F(f1c[1]) isa Nemo.fqPolyRepFieldElem
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]

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

   @test AbstractAlgebra.expressify((x+a*y)^2) isa Expr
end

@testset "Nemo.FqPolyRepFieldElem" begin
   F, _ = Nemo.Native.finite_field(Nemo.next_prime(Nemo.ZZ(10)^50), 2, "a")
   R, _ = polynomial_ring(F, ["x", "y"])
   @test R isa Singular.PolyRing{Singular.n_FieldElem{Nemo.FqPolyRepFieldElem}}

   F, a = Nemo.Native.finite_field(Nemo.ZZ(7), 2, "a")

   R, (x, y) = polynomial_ring(F, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa Singular.n_FieldElem{Nemo.FqPolyRepFieldElem}
   @test F(f1c[1]) isa Nemo.FqPolyRepFieldElem
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]

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

@testset "Nemo.AbsSimpleNumFieldElem" begin
   U, z = Nemo.polynomial_ring(Nemo.QQ, "z")
   K, a = Nemo.number_field(z^3 + 3z + 1, "a")

   R, (x, y) = polynomial_ring(K, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa Singular.n_FieldElem{Nemo.AbsSimpleNumFieldElem}
   @test K(f1c[1]) isa Nemo.AbsSimpleNumFieldElem
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]

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

@testset "Nemo.gfp_fmpz_mod.polynomial_ring" begin

   U = Nemo.Native.GF(Nemo.ZZRingElem(11))

   R, (x, y) = polynomial_ring(U, ["x", "y"])

   wrappedUtype = Singular.n_FieldElem{Singular.FieldElemWrapper{Nemo.FpField, Nemo.FpFieldElem}}

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa wrappedUtype
   @test U(f1c[1]) isa Nemo.FpFieldElem
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]
   @test base_ring(R)(U(3)) isa wrappedUtype
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
   @test inv(f1c[2])*f1c[2] == 1
   @test gcd(f1c[1], f1c[2]) == U(1)
   @test divexact(f1c[2], f1c[1]) == U(3)
   @test f1c[2] - f1c[1] == U(2)
   @test f1c[1] + f1c[2] == U(4)
   @test length(string((x+y)^2)) > 3
   @test hash((x+y)^2) == hash(x^2+2*x*y+y^2)
   @test deepcopy(f1c[1]) == f1c[1]

   I = Ideal(R, x*y+x^3+1, x*y^2+x^2+1)
   @test ngens(std(I)) == 3
end

@testset "Nemo.gfp_fmpz_mod.WeylAlgebra" begin

   U = Nemo.Native.GF(Nemo.ZZRingElem(11))

   R, (x, y, dx, dy) = @inferred WeylAlgebra(U, ["x", "y"])

   wrappedUtype = Singular.n_FieldElem{Singular.FieldElemWrapper{Nemo.FpField, Nemo.FpFieldElem}}

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa wrappedUtype
   @test U(f1c[1]) isa Nemo.FpFieldElem
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]
   @test base_ring(R)(U(3)) isa wrappedUtype
   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2
   @test f1 + U(2) == U(2) + f1
   @test f1 - U(2) == -(U(2) - f1)
   @test U(2)*f1 == f1*U(2)
   @test f1*x == x*f1
   @test f1*dx != dx*f1
   @test deepcopy(f1) == f1
   @test f1*f2 == f2*f1
   @test f1^3*f1^3 == f1^6
   @test inv(f1c[2])*f1c[2] == 1
   @test gcd(f1c[1], f1c[2]) == U(1)
   @test f1c[2] - f1c[1] == U(2)
   @test f1c[1] + f1c[2] == U(4)
   @test length(string((x+y)^2)) > 3
   @test hash((x+y)^2) == hash(x^2+2*x*y+y^2)
   @test deepcopy(f1c[1]) == f1c[1]

   I = Ideal(R, x*y+x^3+1, x*y^2+x^2+1)
   @test ngens(std(I)) == 3
end

@testset "Nemo.zzModRingElem" begin
  R = residue_ring(Nemo.ZZ, 15)[1]
  S, (x, y) = polynomial_ring(R, ["x", "y"])
  p = x*2
  @test string(p) == "2*x"
  @test leading_coefficient(p) == 2
  p = (x*3)*5
  @test string(p) == "0"
  @test length(p) == 0
  @test is_zero(p)

  R = residue_ring(Nemo.ZZ, 7)[1]
  F = Fp(7)
  for i in -10:10
    @test R(one(F)*i) == one(R)*i
    @test F(one(R)*i) == one(F)*i
  end

  S = Nemo.Native.GF(5)
  @test_throws ErrorException S(one(F))
  @test_throws ErrorException F(one(S))
end

@testset "Nemo.fpFieldElem" begin
  R = Nemo.Native.GF(7)
  F = Fp(7)
  for i in -10:10
    @test R(one(F)*i) == one(R)*i
    @test F(one(R)*i) == one(F)*i
  end

  S = Nemo.Native.GF(5)
  @test_throws ErrorException S(one(F))
  @test_throws ErrorException F(one(S))
end

@testset "Nemo.ZZModRingElem" begin
   U = Nemo.residue_ring(Nemo.ZZ, Nemo.ZZRingElem(11))[1]
   R, (x, y) = polynomial_ring(U, ["x", "y"])

   wrappedUtype = Singular.n_RingElem{Singular.RingElemWrapper{Nemo.ZZModRing, Nemo.ZZModRingElem}}

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa wrappedUtype
   @test U(f1c[1]) isa Nemo.ZZModRingElem
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]
   @test base_ring(R)(U(3)) isa wrappedUtype
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
   @test inv(f1c[2])*f1c[2] == 1
   @test gcd(f1c[1], f1c[2]) == U(1)
   @test divexact(f1c[2], f1c[1]) == U(3)
   @test f1c[2] - f1c[1] == U(2)
   @test f1c[1] + f1c[2] == U(4)
   @test length(string((x+y)^2)) > 3
   @test hash((x+y)^2) == hash(x^2+2*x*y+y^2)
   @test deepcopy(f1c[1]) == f1c[1]

   I = Ideal(R, x*y+x^3+1, x*y^2+x^2+1)
   @test ngens(std(I)) == 3
end

@testset "Nemo.NemoField.polynomial_ring" begin
   U = Nemo.Generic.fraction_field(Nemo.ZZ)

   R, (x, y) = polynomial_ring(U, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa Singular.n_FieldElem{AbstractAlgebra.Generic.FracFieldElem{Nemo.ZZRingElem}}
   @test U(f1c[1]) isa AbstractAlgebra.Generic.FracFieldElem{Nemo.ZZRingElem}
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]
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

@testset "Nemo.NemoField.FreeAlgebra" begin
   U = Nemo.Generic.fraction_field(Nemo.ZZ)

   R, (x, y) = @inferred FreeAlgebra(U, ["x", "y"], 13)

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa Singular.n_FieldElem{AbstractAlgebra.Generic.FracFieldElem{Nemo.ZZRingElem}}
   @test U(f1c[1]) isa AbstractAlgebra.Generic.FracFieldElem{Nemo.ZZRingElem}
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]
   @test f1 + 2 == 2 + f1
   @test f1 - 2 == -(2 - f1)
   @test 2*f1 == f1*2
   @test f1 + U(2) == U(2) + f1
   @test f1 - U(2) == -(U(2) - f1)
   @test U(2)*f1 == f1*U(2)
   @test f1*x != x*f1
   @test f1 == (@inferred deepcopy(f1))
   @test f1*f2 != f2*f1
   @test f1^3*f1^3 == f1^6
   @test inv(f1c[2]) == U(1, 3)
   @test gcd(f1c[1], f1c[2]) == U(1)
   @test f1c[2] - f1c[1] == U(2)
   @test f1c[1] + f1c[2] == U(4)
   @test length(string((x+y)^2)) > 3
   @test hash((x+y)^2) == hash(x^2+x*y+y^2+y*x)
   @test deepcopy(f1c[1]) == f1c[1]

   R, (x, y) = FreeAlgebra(Singular.Nemo.QQ, ["x", "y"], 5)
   I = @inferred Ideal(R, x*y+x^3+1, x*y^2+x^2+1)
   @test ngens(std(I)) > 2
end

@testset "Nemo.NemoRing" begin
   U, z = Nemo.polynomial_ring(Nemo.ZZ, "z")

   R, (x, y) = polynomial_ring(U, ["x", "y"])

   f1 = 3x*y + x^2 + 2y
   f2 = y^2 + 1
   f3 = x^2 + 2x + 1

   f1c = [c for c in coefficients(f1)]
   @test f1c[1] isa Singular.n_RingElem{Nemo.ZZPolyRingElem}
   @test U(f1c[1]) isa Nemo.ZZPolyRingElem
   @test !isempty(string(f1c[1]))
   @test leading_coefficient(f1) == f1c[1]

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
