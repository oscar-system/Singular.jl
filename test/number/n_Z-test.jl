function test_n_Z_constructors()
   print("n_Z.constructors...")

   @test elem_type(ZZ) == n_Z
   @test elem_type(Integers) == n_Z
   @test parent_type(n_Z) == Integers
   @test base_ring(ZZ) == Union{}

   typeof(ZZ) <: Nemo.Ring

   a = ZZ()

   @test base_ring(a) == Union{}
   @test parent(a) == ZZ

   @test isa(a, n_Z)

   b = ZZ(123)

   @test isa(b, n_Z)

   c = ZZ(BigInt(123))

   @test isa(c, n_Z)

   d = ZZ(c)

   @test isa(d, n_Z)

   f = ZZ(Nemo.ZZ(123))

   @test isa(f, n_Z)

   println("PASS")
end

function test_n_Z()
   test_n_Z_constructors()

   println("")
end

