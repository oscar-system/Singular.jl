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

function test_n_Z_printing()
   print("n_Z.printing...")

   @test string(ZZ(123)) == "123"

   println("PASS")
end

function test_n_Z_manipulation()
   print("n_Z.manipulation...")

   @test isone(one(ZZ))
   @test iszero(zero(ZZ))
   @test isunit(ZZ(1)) && isunit(ZZ(1))
   @test !isunit(ZZ(2)) && !isunit(ZZ(0)) 

   @test numerator(ZZ(2)) == 2
   @test denominator(ZZ(2)) == 1
   @test denominator(ZZ()) == 1

   @test abs(ZZ(-2)) == 2
   @test abs(ZZ(2)) == 2
   @test abs(ZZ()) == 0

   println("PASS")
end

function test_n_Z_unary_ops()
   print("n_Z.unary_ops...")

   @test -ZZ(123) == -123
   @test -ZZ() == 0

   println("PASS")
end

function test_n_Z_binary_ops()
   print("n_Z.binary_ops...")

   a = ZZ(2)
   b = ZZ(3)

   @test a + b == 5
   @test a - b == -1
   @test a*b == 6

   println("PASS")
end

function test_n_Z_adhoc_binary()
   print("n_Z.adhoc_binary...")

   @test ZZ(2) + 3 == 5
   @test 2 + ZZ(3) == 5
   @test ZZ(2) - 3 == -1
   @test 2 - ZZ(3) == -1
   @test 2*ZZ(3) == 6
   @test ZZ(2)*3 == 6

   println("PASS")
end

function test_n_Z()
   test_n_Z_constructors()
   test_n_Z_printing()
   test_n_Z_manipulation()
   test_n_Z_unary_ops()
   test_n_Z_binary_ops()
   test_n_Z_adhoc_binary()

   println("")
end

