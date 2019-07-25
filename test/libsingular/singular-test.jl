function test_getMPZ()
   print("libSingular.n_GetMPZ...")

   a = Singular.QQ(2//7)
   b = Singular.QQ(1//-3)
   c = Singular.ZZ(-4)

   @test 2 == Singular.libSingular.n_GetMPZ(a.ptr, parent(a).ptr)
   @test -1 == Singular.libSingular.n_GetMPZ(b.ptr, parent(b).ptr)
   @test -4 == Singular.libSingular.n_GetMPZ(c.ptr, parent(c).ptr)

   @test BigInt == typeof(Singular.libSingular.n_GetMPZ(a.ptr, parent(a).ptr))

   println("PASS")
end


function test_singular()
   test_getMPZ()

   println("")
end
