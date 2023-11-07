@testset "PolyRing.degree_bound" begin
   R1, = polynomial_ring(QQ, ["x", ])
   @test degree_bound(R1) >= 0

   R2, = polynomial_ring(QQ, ["x", "y"]; degree_bound = 3)
   @test degree_bound(R2) >= 3

   R3, = polynomial_ring(QQ, [ "x$i" for i in 1:100 ]; degree_bound = 3)
   @test degree_bound(R3) >= 3

   R4, = polynomial_ring(QQ, [ "x$i" for i in 1:12 ]; degree_bound = 5)
   @test degree_bound(R4) >= 5

   R5, = polynomial_ring(QQ, [ "x$i" for i in 1:10 ]; degree_bound = 5)
   @test degree_bound(R5) >= 5

   R6, = polynomial_ring(QQ, [ "x$i" for i in 1:10 ]; degree_bound = 10000000)
   @test degree_bound(R6) >= 10000000
end

@testset "PolyRing.ordering.symbol" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
   I = Ideal(R, x+z, y+z)
   I.isGB = true
   S = fres(I, 0)
   M = S[2]
   @test string(M[1]) == "x*gen(2)-y*gen(1)+z*gen(2)-z*gen(1)"

   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering2 = :comp1max)
   I = Ideal(R, x+z, y+z)
   I.isGB = true
   S = fres(I, 0)
   M = S[2]
   @test string(M[1]) == "x*gen(2)-y*gen(1)-z*gen(1)+z*gen(2)"

   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering = :comp1min,
         ordering2 = :degrevlex)
   I = Ideal(R, x+z, y+z)
   I.isGB = true
   S = fres(I, 0)
   M = S[2]
   @test string(M[1]) == "x*gen(2)+z*gen(2)-y*gen(1)-z*gen(1)"

   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering = :comp1max,
         ordering2 = :degrevlex)
   I = Ideal(R, x+z, y+z)
   I.isGB = true
   S = fres(I, 0)
   M = S[2]
   @test string(M[1]) == "[-y-z,x+z]"
end

@testset "PolyRing.ordering.fancy" begin
   function test_monomials(l)
      for i in 2:length(l)
         @test l[i-1] > l[i]
         @test l[i] < l[rand(1:i-1)]
      end
   end

   function test_ordering(o, s, z, w)
      @test length(o) == length(s)
      j = 0
      for i in o
         j += 1
         @test s[j] == ordering_as_symbol(i)
         @test z[j] == ordering_size(i)
         @test w[j] == ordering_weights(i)
      end
   end

   @test_throws ArgumentError ordering_wp(Int[])
   @test_throws ArgumentError ordering_wp([-1, 2])
   @test_throws ArgumentError ordering_ws(Int[])
   @test_throws ArgumentError ordering_ws([0, 1])
   @test_throws ArgumentError ordering_M([1 2; 3 6])
   @test_throws ArgumentError ordering_M([1 2 3; 3 6 7])

   t = ordering_M([1 2; 3 5]; check = false)*ordering_dp()
   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering = t)
   test_ordering(ordering(R), [:matrix, :degrevlex, :comp1min],
                              [2, 2, 0],
                              [Int[1, 2, 3, 5], Int[], Int[]])
   test_monomials([x1^3*x2^2*x4^3, x1^3*x2^2*x3^2, x1^3*x2^2, x1*x2^3, x1^3*x2,
      x1*x2^2, x1^2*x2*x3, x1^2*x2*x4, x2^2, x1*x2, x2, x1, x3, x4, one(R)])

   t = [1 0 1 0; 0 1 0 -1; 1 0 0 1; 0 1 1 0]
   t = ordering_M(Singular.Nemo.matrix(Singular.Nemo.ZZ, t); check = true)
   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering = t)
   test_ordering(ordering(R), [:matrix, :comp1min],
                              [4, 0],
                              [Int[1,0,1,0, 0,1,0,-1, 1,0,0,1, 0,1,1,0], Int[]])
   test_monomials([x1^2*x2^2*x3^2, x1^2*x2*x3^2, x1^2*x3^2, x1^2*x3^2*x4,
      x1*x2^2*x3^2, x1^2*x2*x3, x1*x3^2, x1*x3^2*x4, x1^2*x2^2, x2^2*x3^2*x4,
      x1*x2*x3, x2*x3^2, x1*x2*x3*x4, x2*x3^2*x4, x1*x3, x3^2, x1*x3*x4,
      x3^2*x4, x1^2*x4^2, x1*x3*x4^2, x3^2*x4^2, x1*x2^2, x2^2*x3, x2^2*x3*x4,
      x1*x2, x2*x3, x1*x2*x4, x2*x3*x4, x1, x3, x1*x2*x4^2, x1*x4, x1*x4^2,
      x2^2, x2^2*x4, x2, x2^2*x4^2, x2*x4, one(R), x2*x4^2, x4, x4^2])


   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering =
                                             ordering_Wp([1, 2])*ordering_ip())
   test_ordering(ordering(R), [:wdeglex, :invlex, :comp1min],
                              [2, 2, 0],
                              [Int[1, 2], Int[], Int[]])
   test_monomials([x2^3, x1*x2^2, x1^2*x2, x2^2*x4, x2^2*x3, x2^2, x1^3,
      x1*x2*x4, x1*x2*x3, x1*x2, x1^2*x4, x1^2*x3, x1^2, x2*x4^2, x2*x3*x4,
      x2*x4, x2*x3^2, x2*x3, x2, x1*x4^2, x1*x3*x4, x1*x4, x1*x3^2, x1*x3,
      x1, x4^3, x3*x4^2, x4^2, x3^2*x4, x3*x4, x4, x3^3, x3^2, x3, one(R)])


   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering =
                                            ordering_is(2)*ordering_wp([1, 2]))
   test_ordering(ordering(R), [:neginvlex, :wdegrevlex, :comp1min],
                              [2, 2, 0],
                              [Int[], Int[1, 2], Int[]])
   test_monomials([x4^3, x3*x4^2, x3^2*x4, x4^2, x3^3, x3*x4, x3^2, x4, x3,
      one(R), x1*x4^2, x1*x3*x4, x1*x3^2, x1*x4, x1*x3, x1, x1^2*x4, x1^2*x3,
      x1^2, x1^3, x2*x4^2, x2*x3*x4, x2*x3^2, x2*x4, x2*x3, x2, x1*x2*x4,
      x1*x2*x3, x1*x2, x1^2*x2, x2^2*x4, x2^2*x3, x2^2, x1*x2^2, x2^3])


   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering =
                                            ordering_ls(2)*ordering_ws([1, 2]))
   test_ordering(ordering(R), [:neglex, :negwdegrevlex, :comp1min],
                              [2, 2, 0],
                              [Int[], Int[1, 2], Int[]])
   test_monomials([one(R), x3, x3^2, x4, x3^3, x3*x4, x3^2*x4, x4^2, x3*x4^2,
      x4^3, x2, x2*x3, x2*x3^2, x2*x4, x2*x3*x4, x2*x4^2, x2^2, x2^2*x3,
      x2^2*x4, x2^3, x1, x1*x3, x1*x3^2, x1*x4, x1*x3*x4, x1*x4^2, x1*x2,
      x1*x2*x3, x1*x2*x4, x1*x2^2, x1^2, x1^2*x3, x1^2*x4, x1^2*x2, x1^3])


   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering =
                                            ordering_Ws([1, 2])*ordering_lp(2))
   test_ordering(ordering(R), [:negwdeglex, :lex, :comp1min],
                              [2, 2, 0],
                              [Int[1, 2], Int[], Int[]])
   test_monomials([x3^3, x3^2*x4, x3^2, x3*x4^2, x3*x4, x3, x4^3, x4^2, x4,
      one(R), x1*x3^2, x1*x3*x4, x1*x3, x1*x4^2, x1*x4, x1, x1^2*x3, x1^2*x4,
      x1^2, x2*x3^2, x2*x3*x4, x2*x3, x2*x4^2, x2*x4, x2, x1^3, x1*x2*x3,
      x1*x2*x4, x1*x2, x1^2*x2, x2^2*x3, x2^2*x4, x2^2, x1*x2^2, x2^3])


   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering =
                                                 ordering_Dp(2)*ordering_Ds(2))
   test_ordering(ordering(R), [:deglex, :negdeglex, :comp1min],
                              [2, 2, 0],
                              [Int[], Int[], Int[]])
   test_monomials([x1^3, x1^2*x2, x1*x2^2, x2^3, x1^2, x1^2*x3, x1^2*x4, x1*x2,
      x1*x2*x3, x1*x2*x4, x2^2, x2^2*x3, x2^2*x4, x1, x1*x3, x1*x4, x1*x3^2,
      x1*x3*x4, x1*x4^2, x2, x2*x3, x2*x4, x2*x3^2, x2*x3*x4, x2*x4^2, one(R),
      x3, x4, x3^2, x3*x4, x4^2, x3^3, x3^2*x4, x3*x4^2, x4^3])


   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering =
                                                                 ordering_ds())
   test_ordering(ordering(R), [:negdegrevlex, :comp1min],
                              [4, 0],
                              [Int[], Int[]])
   test_monomials([one(R), x1, x2, x3, x4, x1^2, x1*x2, x2^2, x1*x3, x2*x3,
      x3^2, x1*x4, x2*x4, x3*x4, x4^2, x1^3, x1^2*x2, x1*x2^2, x2^3, x1^2*x3,
      x1*x2*x3, x2^2*x3, x1*x3^2, x2*x3^2, x3^3, x1^2*x4, x1*x2*x4, x2^2*x4,
      x1*x3*x4, x2*x3*x4, x3^2*x4, x1*x4^2, x2*x4^2, x3*x4^2, x4^3])

   R, (x1, x2, x3, x4) = polynomial_ring(QQ, "x".*string.(1:4), ordering =
                                      ordering_a([1, -1, 1, -1])*ordering_dp())
   test_ordering(ordering(R), [:extraweight, :degrevlex, :comp1min],
                              [4, 4, 0],
                              [Int[1, -1, 1, -1], Int[], Int[]])
   test_monomials([x1^3, x1^2*x3, x1*x3^2, x3^3, x1^3*x2, x1^2*x2*x3,
      x1*x2*x3^2, x1^2, x1*x3, x3^2, x1^3*x2^2, x1^2*x2^2*x3, x1^2*x2,
      x1*x2*x3, x2*x3^2, x1^2*x4, x1*x3*x4, x3^2*x4, x1, x3, x1^3*x2^3,
      x1^2*x2^2, x1*x2^2*x3, x1^2*x2*x4, x1*x2*x3*x4, x1*x2, x2*x3, x1*x4,
      x3*x4, one(R), x1^2*x2^3, x1^2*x2^2*x4, x1*x2^2, x2^2*x3, x1*x2*x4,
      x2*x3*x4, x1*x4^2, x3*x4^2, x2, x4, x1*x2^3, x1*x2^2*x4, x1*x2*x4^2,
      x2^2, x2*x4, x4^2, x2^3, x2^2*x4, x2*x4^2, x4^3])
end
