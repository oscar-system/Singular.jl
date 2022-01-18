@testset "caller.Lib" begin
    R, (x, y) = PolynomialRing(Singular.ZZ, ["x", "y"])
    @test Singular.LibSets.isEqualInt(R, Singular.ZZ(1), Singular.ZZ(1)) == 1

    R, (x, y) = PolynomialRing(Singular.QQ, ["x", "y"])

    # Input tests
    @test 1 == Singular.LibSets.isEqualInt(R, 1, 1)
    @test 1 == Singular.LibSets.isEqualInt(R, x, x)
    @test 1 == Singular.LibSets.isEqualInt(R, "aa", "aa")
    @test 1 == Singular.LibSets.isEqualInt(R, Singular.QQ(1), Singular.QQ(1))
    @test 1 == Singular.LibSets.isEqualInt(R, [1,2,3], [1,2,3])
    @test 1 == Singular.LibSets.isEqualInt(R, [1 2; 3 4], [1 2; 3 4])
    @test 1 == Singular.LibSets.isEqualInt(R, BigInt[1 2; 3 4],
                                              Nemo.matrix(Nemo.ZZ, [1 2; 3 4]))

    R, (x,y,z) = PolynomialRing(Singular.QQ, ["x", "y", "z"])

    i1 = Singular.LibPolylib.cyclic(R, 3)
    i2 = Ideal( R, x+y+z, x*y+x*z+y*z, x*y*z-1 )
    @test equal(i1, i2)

    vec = FreeModule(R,2)([x,y])
    mod = Singular.Module(R, vec)
    i1 = Singular.LibPolylib.mod2id(R,mod,[1,2])
    i2 = Ideal(R, x^2, x*y, y^2, x^2 )
    @test equal(i1, i2)

    i1 = Ideal(R, x, y)
    i2 = Ideal(R, x^2, x*y, y^2, x, y)
    mod = Singular.LibPolylib.id2mod(R, i1, [1,2])
    i1 = Singular.LibPolylib.mod2id(R, mod, [1,2])
    @test equal(i1,i2)
    @test Singular.LibPolylib.content(R,vec) == 1
    @test Singular.LibPolylib.lcm(R,x) == x

    i1 = Ideal(R, x*z, y*z, x^3-y^3)
    @test Singular.LibStandard.res(R,i1,0) isa Singular.sresolution
    i1 = Ideal(R, x*z, y*z, x^3-y^3)
    @test Singular.LibPrimdec.primdecGTZ(R,i1) isa Array

    i1 = Ideal(R, x*z, y*z, x^3-y^3)
    @test Singular.LibNormal.normal(R, i1, "withDelta", "prim") isa Array

    @test Singular.LibNormal.normal(i1, "withDelta", "prim") isa Array

    R, (x, y, z) = PolynomialRing(Singular.QQ, ["x", "y", "z"])
    m = Singular.Matrix(R, [0 1 2; 2 1 0])
    id = Singular.LibSing4ti2.markov4ti2(m)
    id = std(id)
    @test 0 == reduce(x*z-y^2, id)

    R, (x, y, z, w) = PolynomialRing(Singular.QQ, ["x", "y", "z", "w"])
    m = Singular.Matrix(R, [0 1 2 3; 3 2 1 0])
    id = Singular.LibSing4ti2.graver4ti2(m)
    id = std(id)
    @test 0 == reduce(x*z-y^2, id)
    @test 0 == reduce(x^2*w-y^3, id)
    @test 0 == reduce(x*w-y*z, id)
    @test 0 == reduce(y*w-z^2, id)
    @test 0 == reduce(x*w^2-z^3, id)

    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9) =
                      PolynomialRing(Singular.QQ, ["x"*string(i) for i in 1:9])
    m = Singular.Matrix(R, [1  1  1 -1 -1 -1  0  0  0;
                            1  1  1  0  0  0 -1 -1 -1;
                            0  1  1 -1  0  0 -1  0  0;
                            1  0  1  0 -1  0  0 -1  0;
                            1  1  0  0  0 -1  0  0 -1;
                            0  1  1  0 -1  0  0  0 -1;
                            1  1  0  0 -1  0 -1  0  0])
    id = Singular.LibSing4ti2.hilbert4ti2(m)
    id = std(id)
    @test 0 == reduce(x1^2*x3*x5*x6^2*x7*x8^2 - 1, id)
    @test 0 == reduce(x1*x3^2*x4^2*x5*x8^2*x9 - 1, id)
    @test 0 == reduce(x2^2*x3*x4^2*x5*x7*x9^2 - 1, id)
    @test 0 == reduce(x1*x2^2*x5*x6^2*x7^2*x9 - 1, id)
    @test 0 == reduce(x1*x2*x3*x4*x5*x6*x7*x8*x9 - 1, id)
end

@testset "caller.LibFinvar" begin
   Qw, (w,) = FunctionField(QQ, ["w"])
   K, w = AlgebraicExtensionField(Qw, w^2+w+1)
   R, (x, y, z) = PolynomialRing(K, ["x","y","z"])
   M1 = Matrix(R, [0 1 0; 0 0 1; 1 0 0])
   M2 = Matrix(R, [1 0 0; 0 w 0; 0 0 w^2])
   P, S, IS = Singular.LibFinvar.invariant_ring(M1, M2)
   @test P isa Singular.smatrix{Singular.spoly{Singular.n_algExt}}
   @test typeof(P) == typeof(S) == typeof(IS)
   @test nrows(P) == 1
   @test ncols(P) == 3
end

@testset "caller.LibSolve" begin
   R, (x,) = PolynomialRing(Singular.QQ, ["s"])
   l = Singular.LibSolve.solve(Ideal(R, [x^5+1]), "nodisplay")
   S = l[1]
   @test S isa Singular.PolyRing{Singular.n_unknownsingularcoefficient}
   for y in l[2][:SOL]
      @test isa(y^5+1, Singular.n_unknownsingularcoefficient)
   end
end

@testset "caller.LibNormal" begin
   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
   I = Ideal(R, z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test S isa Singular.PolyRing{n_Q}
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]) == QQ

   F, a = FiniteField(5, 3, "b")
   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test S isa Singular.PolyRing{n_GF}
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == F

   R, (x, y, z) = PolynomialRing(Fp(17), ["x", "y", "z"])
   I = Ideal(R, z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test S isa Singular.PolyRing{n_Zp}
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]) == Fp(17)

   F, (a, b, c) = FunctionField(Fp(17), ["a", "b", "c"])
   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, b*c*z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test S isa Singular.PolyRing{n_transExt}
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == F

   F, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])
   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, b*c*z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test S isa Singular.PolyRing{n_transExt}
   @test base_ring(gens(l[1][1][2][:norid])[1]) == F
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == F

   F, (Fa,) = FunctionField(QQ, ["a"])
   K, a = AlgebraicExtensionField(F, Fa^2 + 1)
   R, (x, y, z) = PolynomialRing(K, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test S isa Singular.PolyRing{n_algExt}
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == K

   F, (Fa,) = FunctionField(Fp(7), ["a"])
   K, a = AlgebraicExtensionField(F, Fa^2 + 1)
   R, (x, y, z) = PolynomialRing(K, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test S isa Singular.PolyRing{n_algExt}
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == K
end

@testset "caller.lists" begin
    r, (x, y) = PolynomialRing(QQ, ["x", "y"])
    f = y^2*(y-1)^3-x^5
    A1 = Any[Any[Ideal(r, [x,y-1]), 1],2]
    A2 = Any[Any[Ideal(r, [y^2-y+1, x]), 1],3]
    A3 = Any[Any[Ideal(r, [r(1), y-x]), x],-2]
    D = Any[A1, A2]
    E = Any[A3]
    L = Singular.LibHess.RiemannRochHess(r, f, Any[D, E], "free")

    @test length(L[1]) == 5
    @test length(findall(P->P== 5*x^2*y^3-5*x*y^4-5*x^2*y^2+10*x*y^3-2*y^4-5*x*y^2+4*y^3-2*y^2, L[1])) == 1
    @test length(findall(P->P== 2*x^3*y^2+x^2*y^3-3*x*y^4-2*x^3*y-x^2*y^2+6*x*y^3-3*x*y^2, L[1])) == 1
    @test length(findall(P->P== 5*x^3*y^2-5*x*y^4-4*x^3*y+10*x*y^3-5*x*y^2, L[1])) == 1
    @test length(findall(P->P== 5*x^4*y-5*x*y^4-x^3*y+10*x*y^3-5*x*y^2, L[1])) == 1
    @test length(findall(P->P== x^5+2*x^4*y-3*x*y^4+6*x*y^3-3*x*y^2, L[1])) == 1
    @test L[2] == x^5

    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
    F = y^2*(y-z)-x^3
    I = Ideal(R, [x, y-z])
    J = Ideal(R, [x-z, y])
    M = Singular.LibHess.RiemannRochHess(R, F, Any[I, J], "ideals")

    @test length(M[1]) == 2
    @test length(findall(P->P== y, M[1])) == 1
    @test length(findall(P->P== x, M[1])) == 1
    @test M[2] == x
end

@testset "caller.lookup_library_symbol" begin
    F = FiniteField(3, 1, "a")[1]
    R, (x, y, z) = PolynomialRing(F, ["x", "y", "z" ])
    A = Matrix(R, [0 1 0; -1 0 0; 0 0 -1])
    Singular.LibFinvar.reynolds_molien(A, "")

    x = Singular.lookup_library_symbol("Finvar", "newring")
    @test length(x) == 2
    @test x[1] isa PolyRing
    @test x[2][:M] isa smatrix

    @test_throws Exception Singular.lookup_library_symbol("Finvar", "meiyou")
    @test_throws Exception Singular.lookup_library_symbol("blah", "bluh")
end

@testset "caller.noncommutative" begin
    F, (q1, q2) = FunctionField(QQ, ["q1", "q2"])
    R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])
    S, (x, y, z) = GAlgebra(R, Singular.Matrix(R, [0 q2 q1; 0 0 1; 0 0 0]),
                               Singular.Matrix(R, [0 x z; 0 0 0; 0 0 0]))
    @test Singular.LibNctools.ndcond(S) isa sideal

    AA, (x, y, z, t) = PolynomialRing(QQ, ["x", "y", "z", "t"])
    D = zero_matrix(AA, 4, 4)
    D[1,2] = -z; D[1,3] = 2*x; D[2,3] = -2*y
    A, (x, y, z, t) = GAlgebra(AA, Singular.Matrix(AA, [0 1 1 1; 0 0 1 1; 0 0 0 1; 0 0 0 1]), D)
    @test Singular.LibCentral.center(A, 3) isa sideal
end
