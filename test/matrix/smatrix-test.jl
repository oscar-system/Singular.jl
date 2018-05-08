function test_smatrix_constructors()
   print("smatrix.constructors...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x, y)
   M = Singular.Matrix(I)
   S = parent(M)

   @test elem_type(S) == smatrix{spoly{n_Q}}
   @test elem_type(MatrixSpace{spoly{n_Q}}) == smatrix{spoly{n_Q}}
   @test parent_type(smatrix{spoly{n_Q}}) == MatrixSpace{spoly{n_Q}}

   @test base_ring(S) == R
   @test base_ring(M) == R

   @test typeof(S) <: AbstractAlgebra.Set
   @test typeof(M) <: AbstractAlgebra.SetElem

   N = std(I)
   V = Singular.Matrix(N)

   @test isa(M, smatrix)
   @test isa(V, smatrix)

   println("PASS")
end

function test_smatrix_manipulation()
   print("smatrix.manipulation...")

   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x, y)
   M = Singular.Matrix(I)

   @test nrows(M) == 1
   @test ncols(M) == 2

   @test M[1, 1] == x
   @test M[1, 2] == y

   @test !iszero(M)
   @test iszero(M - M)

   println("PASS")
end

function test_smatrix()
   test_smatrix_constructors()
   test_smatrix_manipulation()

   println("")
end

