@testset "caller.Lib..." begin
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

    # for 4ti2 tests:
    # m = matrix(R, [0 1 2; 2 1 0]) does not work
    # furthermore, matrices are not implemented properly in caller.jl
    # so we pass in a module and let singular do the conversion to a matrix

    R, (x, y, z) = PolynomialRing(Singular.QQ, ["x", "y", "z"])
    m = Singular.smodule{n_Q}(R, vector(R, R(0), R(2)),
                                 vector(R, R(1), R(1)),
                                 vector(R, R(2), R(0)))
    id = Singular.LibSing4ti2.markov4ti2(m)
    id = std(id)
    @test 0 == reduce(x*z-y^2, id)

    R, (x, y, z, w) = PolynomialRing(Singular.QQ, ["x", "y", "z", "w"])
    m = Singular.smodule{n_Q}(R, vector(R, R(0), R(3)),
                                 vector(R, R(1), R(2)),
                                 vector(R, R(2), R(1)),
                                 vector(R, R(3), R(0)))
    id = Singular.LibSing4ti2.graver4ti2(m)
    id = std(id)
    @test 0 == reduce(x*z-y^2, id)
    @test 0 == reduce(x^2*w-y^3, id)
    @test 0 == reduce(x*w-y*z, id)
    @test 0 == reduce(y*w-z^2, id)
    @test 0 == reduce(x*w^2-z^3, id)

    R, (x1, x2, x3, x4, x5, x6, x7, x8, x9) =
                      PolynomialRing(Singular.QQ, ["x"*string(i) for i in 1:9])
    m = Singular.smodule{n_Q}(R,
                        vector(R, R(1), R(1), R(0), R(1), R(1), R(0), R(1)),
                        vector(R, R(1), R(1), R(1), R(0), R(1), R(1), R(1)),
                        vector(R, R(1), R(1), R(1), R(1), R(0), R(1), R(0)),
                        vector(R, R(-1),R(0), R(-1),R(0), R(0), R(0), R(0)),
                        vector(R, R(-1),R(0), R(0), R(-1),R(0), R(-1),R(-1)),
                        vector(R, R(-1),R(0), R(0), R(0), R(-1),R(0), R(0)),
                        vector(R, R(0), R(-1),R(-1),R(0), R(0), R(0), R(-1)),
                        vector(R, R(0), R(-1),R(0), R(-1),R(0), R(0), R(0)),
                        vector(R, R(0), R(-1),R(0), R(0), R(-1),R(-1),R(0)))
    id = Singular.LibSing4ti2.hilbert4ti2(m)
    id = std(id)
    @test 0 == reduce(x1^2*x3*x5*x6^2*x7*x8^2 - 1, id)
    @test 0 == reduce(x1*x3^2*x4^2*x5*x8^2*x9 - 1, id)
    @test 0 == reduce(x2^2*x3*x4^2*x5*x7*x9^2 - 1, id)
    @test 0 == reduce(x1*x2^2*x5*x6^2*x7^2*x9 - 1, id)
    @test 0 == reduce(x1*x2*x3*x4*x5*x6*x7*x8*x9 - 1, id)
end

@testset "caller.LibSolve" begin
   R, (x,) = PolynomialRing(Singular.QQ, ["s"])
   l = Singular.LibSolve.solve(Ideal(R, [x^5+1]))
   S = l[1]
   @test isa(S, Singular.PolyRing{Singular.n_unknownsingularcoefficient})
   for y in l[2][:SOL]
      @test isa(y^5+1, Singular.n_unknownsingularcoefficient)
   end
end

@testset "caller.LibNormal..." begin
   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
   I = Ideal(R, z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test isa(S, Singular.PolyRing{n_Q})
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]) == QQ

   F, a = FiniteField(5, 3, "b")
   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test isa(S, Singular.PolyRing{n_GF})
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == F

   R, (x, y, z) = PolynomialRing(Fp(17), ["x", "y", "z"])
   I = Ideal(R, z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test isa(S, Singular.PolyRing{n_Zp})
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]) == Fp(17)

   F, (a, b, c) = FunctionField(Fp(17), ["a", "b", "c"])
   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, b*c*z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test isa(S, Singular.PolyRing{n_transExt})
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == F

   F, (a, b, c) = FunctionField(QQ, ["a", "b", "c"])
   R, (x, y, z) = PolynomialRing(F, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, b*c*z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   @test isa(S, Singular.PolyRing{n_transExt})
   @test base_ring(gens(l[1][1][2][:norid])[1]) == F
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == F

   F, (Fa,) = FunctionField(QQ, ["a"])
   K, a = AlgebraicExtensionField(F, Fa^2 + 1)
   R, (x, y, z) = PolynomialRing(K, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   # the correct tests for when the possible bug is fixed
   #@test isa(S, Singular.PolyRing{n_algExt})
   #@test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == K
   @test isa(S, Singular.PolyRing{n_transExt})
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*Fa) == F

   F, (Fa,) = FunctionField(Fp(7), ["a"])
   K, a = AlgebraicExtensionField(F, Fa^2 + 1)
   R, (x, y, z) = PolynomialRing(K, ["x", "y", "z"])
   I = Ideal(R, a*z-x^4, z-y^6)
   l = Singular.LibNormal.normal(I)
   S = l[1][1][1]
   # the correct tests for when the possible bug is fixed
   #@test isa(S, Singular.PolyRing{n_algExt})
   #@test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*a) == K
   @test isa(S, Singular.PolyRing{n_transExt})
   @test base_ring(gens(l[1][1][2][:norid])[1]*gens(S)[1]*Fa) == F

@testset "caller.lists..." begin
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
