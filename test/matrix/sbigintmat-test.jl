@testset "sbigintmat.constructors" begin
   m = Singular.sbigintmat(2, 3)
   @test size(m) == (2, 3)
   @test nrows(m) == 2
   @test ncols(m) == 3
   @test length(string(m)) > 6
end

@testset "sbigintmat.conversion" begin
   a = Nemo.matrix(Nemo.ZZ, [1 2 3; 4 5 6])
   @test a == Nemo.matrix(Nemo.ZZ, Singular.sbigintmat(a))

   a = BigInt[1 2 3; 4 5 6]
   @test a == Array(Singular.sbigintmat(a))

   R = Nemo.matrix_ring(Nemo.ZZ, 2)
   a = R([1 2; 3 4])
   b = R([5 6; 7 8])
   @test a == R(Singular.sbigintmat(a))
   @test a*b == R(Singular.sbigintmat(a*b))
end
