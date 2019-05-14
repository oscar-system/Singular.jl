function test_caller()
    print("caller...")

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
    
    i1 = Singular.LibPoly.cyclic(R, 3)
    i2 = Ideal( R, x+y+z, x*y+x*z+y*z, x*y*z-1 )
    @test equal(i1, i2)

    vec = FreeModule(R,2)([x,y])
    mod = Singular.Module(R, vec)
    i1 = Singular.LibPoly.mod2id(R,mod,[1,2])
    i2 = Ideal(R, x^2, x*y, y^2, x^2 )
    @test equal(i1, i2)

    i1 = Ideal(R, x, y)
    i2 = Ideal(R, x^2, x*y, y^2, x, y)
    mod = Singular.LibPoly.id2mod(R, i1, [1,2])
    i1 = Singular.LibPoly.mod2id(R, mod, [1,2])
    @test equal(i1,i2)
    @test Singular.LibPoly.content(R,vec) == 1
    @test Singular.LibPoly.lcm(R,x) == x

    i1 = Ideal(R, x*z, y*z, x^3-y^3)
    @test Singular.LibStandard.res(R,i1,0) isa Singular.sresolution
    i1 = Ideal(R, x*z, y*z, x^3-y^3)
    @test Singular.LibPrimdec.primdecGTZ(R,i1) isa Array

    i1 = Ideal(R, x*z, y*z, x^3-y^3)
    @test Singular.LibNormal.normal(i1, "withDelta", "prim") isa Array

    println("PASS")
end
