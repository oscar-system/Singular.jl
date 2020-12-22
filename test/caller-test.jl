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
