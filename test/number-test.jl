testrings = [
      (n_Z, Integers, Union{}, ZZ),
      (n_Q, Rationals, ZZ, QQ),
      (n_Zn, N_ZnRing, ZZ, ResidueRing(ZZ, 7)),
      (n_Zp, N_ZpField, Union{}, Fp(7)),
      (n_GF, N_GField, Union{}, FiniteField(7, 2, "x")[1]),
      (n_transExt, N_FField, QQ, FunctionField(QQ, ["a", "b", "c"])[1]),
      ]


@testset "number.parent and type relations for $(elemT)..." for (elemT, ringT, parentRing, R) in testrings
   @test elem_type(R) == elemT
   @test elem_type(ringT) == elemT
   @test parent_type(elemT) == ringT
   @test base_ring(R) == parentRing
   @test R isa ringT
end

@testset "number.integer_constructors for $(elemT)..." for (elemT, ringT, parentRing, R) in testrings
   a = R()
   @test base_ring(a) == parentRing
   @test parent(a) == R
   @test a isa elemT

   # "small" integers
   b = R(3)
   @test b isa elemT

   @testset "small integer for $(elemT) and $U..." for U in (BigInt, ZZ, R, Nemo.ZZ)
      c = R(U(3))
      @test c isa elemT
      @test b == c
   end

   # "big" integers
   b = R(3)^80
   @test b isa elemT

   @testset "big integer for $(elemT) and $U..." for U in (BigInt, ZZ, R, Nemo.ZZ)
      c = R(U(3)^80)
      @test c isa elemT
      @test b == c
   end
end

# test adhoc arithmetic resp. arithmetic via promotion_rules
@testset "number.adhoc_binary for $(elemT)..." for (elemT, ringT, parentRing, R) in testrings
   @testset "number.adhoc_binary for $(elemT) and $T..." for T in (Int8, UInt8, Int, BigInt, Nemo.ZZ, ZZ)

      a = R(2)
      b = T(3)

      @test a + b == R(5)
      @test a - b == R(-1)
      @test a * b == R(6)

      a = T(2)
      b = R(3)

      @test a + b == R(5)
      @test a - b == R(-1)
      @test a * b == R(6)
   end
end

include("number/n_Z-test.jl")
include("number/n_Zp-test.jl")
include("number/n_Zn-test.jl")
include("number/n_GF-test.jl")
include("number/n_transExt-test.jl")
include("number/n_Q-test.jl")
