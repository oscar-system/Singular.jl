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

   @test S isa AbstractAlgebra.Set
   @test M isa AbstractAlgebra.SetElem

   @test S(M) === M

   N = std(I)
   V = Singular.Matrix(N)

   Z = zero_matrix(R, 2 ,3)
   Z2 = S()
   E = identity_matrix(R, 3)

   @test isa(E, smatrix)
   @test isa(M, smatrix)
   @test isa(V, smatrix)
   @test isa(Z, smatrix)
   @test iszero(Z)
   @test isa(Z2, smatrix)
   @test iszero(Z2)

   D = S(3) # TODO: test more interesting matrix dimensions
   @test D isa smatrix
   @test D[1, 1] == 3
   @test D[1, 2] == 0

   D = S(x+y)
   @test D isa smatrix
   @test D[1, 1] == x+y
   @test D[1, 2] == 0

   Rx, _ = PolynomialRing(QQ, ["x"])
   Z = zero_matrix(Rx, 2, 3)
   @test_throws ArgumentError S(Z)
   @test_throws Exception S(Z[1, 1])
end

@testset "smatrix.manipulation..." begin
   R, (x, y) = PolynomialRing(QQ, ["x", "y"])

   I = Ideal(R, x, y)
   M = Singular.Matrix(I)

   Z0 = zero_matrix(R, 3, 3)
   Z1 = zero_matrix(R, 3, 2)
   Z2 = zero_matrix(R, 2, 3)

   Z0[1, 1] = x

   @test nrows(M) == 1
   @test ncols(M) == 2

   @test M[1, 1] == x
   @test M[1, 2] == y
   @test Z0[1, 1] == x

   @test !iszero(M)
   @test iszero(M - M)
   @test transpose(Z1) == Z2
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
