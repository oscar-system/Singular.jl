@testset "PolyRing.degree_bound..." begin
   R1, = PolynomialRing(QQ, ["x", ])
   @test degree_bound(R1) == 65535

   R2, = PolynomialRing(QQ, ["x", "y"]; degree_bound = 3)
   @test degree_bound(R2) == 65535

   R3, = PolynomialRing(QQ, [ "x$i" for i in 1:100 ]; degree_bound = 3)
   @test degree_bound(R3) == 3

   R4, = PolynomialRing(QQ, [ "x$i" for i in 1:12 ]; degree_bound = 5)
   @test degree_bound(R4) == 31

   R5, = PolynomialRing(QQ, [ "x$i" for i in 1:10 ]; degree_bound = 5)
   if Sys.WORD_SIZE == 64
      @test degree_bound(R5) == 63
   elseif Sys.WORD_SIZE == 32
      @test degree_bound(R5) == 7
   end

   R6, = PolynomialRing(QQ, [ "x$i" for i in 1:10 ]; degree_bound = 10000000)
   if Sys.WORD_SIZE == 64
      @test degree_bound(R6) == 4294967295
   elseif Sys.WORD_SIZE == 32
      @test degree_bound(R6) == 2147483647
   end
end

@testset "PolyRing.ordering..." begin
   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
   I = Ideal(R, x+z, y+z)
   I.isGB = true
   S = fres(I, 0)
   M = S[2]
   @test string(M[1]) == "x*gen(2)-y*gen(1)+z*gen(2)-z*gen(1)"

   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering2 = :comp1max)
   I = Ideal(R, x+z, y+z)
   I.isGB = true
   S = fres(I, 0)
   M = S[2]
   @test string(M[1]) == "x*gen(2)-y*gen(1)-z*gen(1)+z*gen(2)"

   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering = :comp1min,
         ordering2 = :degrevlex)
   I = Ideal(R, x+z, y+z)
   I.isGB = true
   S = fres(I, 0)
   M = S[2]
   @test string(M[1]) == "x*gen(2)+z*gen(2)-y*gen(1)-z*gen(1)"

   R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering = :comp1max,
         ordering2 = :degrevlex)
   I = Ideal(R, x+z, y+z)
   I.isGB = true
   S = fres(I, 0)
   M = S[2]
   @test string(M[1]) == "[-y-z,x+z]"
end
