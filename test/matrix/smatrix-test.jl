@testset "smatrix.constructors..." begin
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
end

@testset "smatrix.manipulation..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x, y)
   M = Singular.Matrix(I)

   @test nrows(M) == 1
   @test ncols(M) == 2

   @test M[1, 1] == x
   @test M[1, 2] == y

   @test !iszero(M)
   @test iszero(M - M)

   @test deepcopy(M) == M
end

@testset "smatrix.binary_ops..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I1 = Ideal(R, x, y)
   I2 = Ideal(R, x*y + 1, x^2 + 1)
   I3 = Ideal(R, x^2*y^2 + x + 1)

   M1 = Singular.Matrix(I1)
   M2 = Singular.Matrix(I2)
   M3 = Singular.Matrix(I3)

   @test M1 + M2 == M2 + M1
   @test (M1 + M2) - M2 == M1
   @test M3*(M1 + M2) == M3*M1 + M3*M2
end

@testset "smatrix.comparison..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I1 = Ideal(R, x, y)
   I2 = Ideal(R, x*y + 1, x^2 + 1)

   M1 = Singular.Matrix(I1)
   M2 = Singular.Matrix(I2)

   @test M1 != M2
   @test M1 == M1
   @test M1 == deepcopy(M1)
end

