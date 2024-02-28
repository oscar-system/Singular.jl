@testset "sideal.constructors" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Ideal(R, x, y)
   S = parent(I)

   @test elem_type(S) == sideal{spoly{n_Q}}
   @test elem_type(IdealSet{spoly{n_Q}}) == sideal{spoly{n_Q}}
   @test parent_type(sideal{spoly{n_Q}}) == IdealSet{spoly{n_Q}}

   @test base_ring(S) == R
   @test base_ring(I) == R

   @test S isa AbstractAlgebra.Set
   @test I isa sideal

   I1 = @inferred Ideal(R)
   I2 = @inferred Ideal(R, x + y)
   I3 = @inferred Ideal(R, x + y, y^2 + 2)
   I4 = @inferred Ideal(R, [x + y, y^2 + 2])
   I5 = @inferred MaximalIdeal(R, 0)
   I6 = @inferred MaximalIdeal(R, 4)

   @test isa(I1, sideal)
   @test isa(I2, sideal)
   @test isa(I3, sideal)
   @test isa(I4, sideal)
   @test isa(I5, sideal)
   @test isa(I6, sideal)

   @test_throws DomainError MaximalIdeal(R, -rand(1:99))
   if sizeof(Cint) < sizeof(Int)
      @test_throws DomainError MaximalIdeal(R, typemax(Int))
   end
end

@testset "sideal.spoly.manipulation" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I0 = Ideal(R)
   I1 = Ideal(R, x, y)

   @test ngens(I0) == 1
   @test ngens(I1) == 2
   @test x in gens(I1) && y in gens(I1)

   I2 = @inferred deepcopy(I1)

   @test isequal(I1, I2)
   normalize!(I1)
   @test isequal(I1, I2)

   @test I2[2] == y
   I2[2] = x^2 + 1
   @test I2[2] == x^2 + 1

   @test iszero(I0)
   @test is_zerodim(std(I1))
   @test dimension(std(I0)) == 2
   @test dimension(std(I1)) == 0

   @test is_constant(Ideal(R, R(1), R(2)))

   @test is_var_generated(Ideal(R, x))
   @test is_var_generated(Ideal(R, y))
   @test is_var_generated(Ideal(R, x, y))
   @test !is_var_generated(Ideal(R, R(1)))
   @test !is_var_generated(Ideal(R, x + y))

   R, (x, y) = polynomial_ring(ZZ, ["x", "y"])
   @test -1 == dimension(std(Ideal(R, R(-1))))
   @test -1 == dimension(std(Ideal(R, x, R(1))))
   @test 0 == dimension(std(Ideal(R, x, y, R(3))))
   @test 1 == dimension(std(Ideal(R, x, R(3))))
   @test 1 == dimension(std(Ideal(R, x - y, R(2))))
   @test 2 == dimension(std(Ideal(R, R(5))))
   @test 2 == dimension(std(Ideal(R, R(35))))

   R, (x, y, z) = polynomial_ring(Fp(32003), ["x", "y", "z"], ordering = ordering_ds())
   a, b, c, t = 11, 10, 3, 1
   f = x^a+y^b+z^(3*c)+x^(c+2)*y^(c-1)+x^(c-1)*y^(c-1)*z^3+x^(c-2)*y^c*(y^2+t*x)^2
   i = @inferred jacobian_ideal(f)
   i0 = @inferred std(i)
   @test degree(i0) == (0, 314)

   R, (x, y) = polynomial_ring(Fp(32003), ["x", "y"], ordering = ordering_ds())
   f = (x^3+y^5)^2+x^2*y^7
   @test mult(std(@inferred jacobian_ideal(f))) == 46
   @test mult(std(Ideal(R, f))) == 6

   R, (x, y) = polynomial_ring(QQ, ["x", "y"]; ordering= :lex)
   @test_throws ErrorException iszerodim(Ideal(R, x+y,x^2))
end

@testset "sideal.spoly.homogeneous" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
   I = Ideal(R, x^2+y^3+z^4, x^5+y^4+z^3)
   @test !is_homogeneous(I)
   @test is_homogeneous(homogenize(I, z))
   I = Ideal(R, x^2+y^2+z^2, x^3+y^3+z^3)
   @test is_homogeneous(I)
   @test is_homogeneous(homogenize(I, z))

   R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], ordering = ordering_wp([6,4,3,1]))
   I = Ideal(R, x^2+y^3+z^4, x^1+z^2)
   @test is_homogeneous(I)
   @test is_homogeneous(homogenize(I, w))
   I = Ideal(R, x^2+y^2+z^2, x^3+y^3+z^3)
   @test !is_homogeneous(I)
   @test is_homogeneous(homogenize(I, w))
end

@testset "sideal.spluralg.manipulation" begin

   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x + y; 0 0]))

   I0 = Ideal(R)
   I1 = Ideal(R, x, y)

   @test ngens(I0) == 1
   @test ngens(I1) == 2
   @test x in gens(I1) && y in gens(I1)

   I2 = @inferred deepcopy(I1)

   @test isequal(I1, I2)
   normalize!(I1)
   @test isequal(I1, I2)

   @test I2[2] == y
   I2[2] = x^2 + 1
   @test I2[2] == x^2 + 1

   @test iszero(I0)

   @test is_constant(Ideal(R, R(1), R(2)))

   @test is_var_generated(Ideal(R, x))
   @test is_var_generated(Ideal(R, y))
   @test is_var_generated(Ideal(R, x, y))
   @test !is_var_generated(Ideal(R, R(1)))
   @test !is_var_generated(Ideal(R, x + y))

end

@testset "sideal.slpalg.manipulation" begin

   R, (x, y) = FreeAlgebra(Fp(11), ["x", "y"], 9)

   I0 = Ideal(R)
   I1 = Ideal(R, x, y)

   @test ngens(I0) == 1
   @test ngens(I1) == 2
   @test x in gens(I1) && y in gens(I1)

   I2 = @inferred deepcopy(I1)

   @test isequal(I1, I2)
   normalize!(I1)
   @test isequal(I1, I2)

   @test I2[2] == y
   I2[2] = x^2 + 1
   @test I2[2] == x^2 + 1

   @test iszero(I0)

   @test is_constant(Ideal(R, R(1), R(2)))

   @test is_var_generated(Ideal(R, x))
   @test is_var_generated(Ideal(R, y))
   @test is_var_generated(Ideal(R, x, y))
   @test !is_var_generated(Ideal(R, R(1)))
   @test !is_var_generated(Ideal(R, x + y))

end

@testset "sideal.binary_ops" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I1 = @inferred Ideal(R, x)
   I2 = @inferred Ideal(R, y)
   I3 = @inferred Ideal(R, x, y)
   I4 = @inferred Ideal(R, x*y)

   @test equal(I1 + I2, I3)
   @test equal(I1*I2, I4)

   @test equal(I3*x, Ideal(R, x^2, y*x))
   @test equal(y*I3, Ideal(R, y^2, y*x))
   @test equal(I3*2, I3)
   @test equal(2*I3, I3)

   R, (x, dx) = WeylAlgebra(QQ, ["x"])

   I1 = @inferred Ideal(R, x)
   I2 = @inferred Ideal(R, dx)

   @test equal(I1*x, x*I1)
   @test !equal(I1*dx, dx*I1)
   @test !equal(I2*x, x*I2)
   @test equal(I2*dx, dx*I2)
end

@testset "sideal.powering" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Ideal(R, x^2, x*y + 1)

   @test equal(I^0, Ideal(R, R(1)))

   S = I

   for i = 1:5
      @test equal(S, I^i)
      S *= I
   end

   @test_throws DomainError I^(-rand(1:99))
   if sizeof(Cint) < sizeof(Int)
      @test_throws DomainError I^typemax(Int)
   end
end

@testset "sideal.contains" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   @test contains(Ideal(R, x, y), Ideal(R, x))
   @test !contains(Ideal(R, x), Ideal(R, x, y))

   r, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x + y; 0 0]))

   @test contains(Ideal(R, x, y), Ideal(R, x))
   @test !contains(Ideal(R, x), Ideal(R, x, y))

   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 10)

   @test contains(Ideal(R, x, y), Ideal(R, x))
   @test !contains(Ideal(R, x), Ideal(R, x, y))
end

@testset "sideal.comparison" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   @test equal(Ideal(R, x, y, x - y), Ideal(R, x, y, x + y))
   @test !equal(Ideal(R, x), Ideal(R, x, y))
   @test !isequal(Ideal(R, x, y, x - y), Ideal(R, x, y, x + y))
   @test isequal(Ideal(R, x, y), Ideal(R, x, y))

   r, (x, y) = polynomial_ring(Fp(7), ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x + y; 0 0]))

   @test equal(Ideal(R, x, y, x - y), Ideal(R, x, y, x + y))
   @test !equal(Ideal(R, x), Ideal(R, x, y))
   @test !isequal(Ideal(R, x, y, x - y), Ideal(R, x, y, x + y))
   @test isequal(Ideal(R, x, y), Ideal(R, x, y))

   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 10)

   @test equal(Ideal(R, x, y, x - y), Ideal(R, x, y, x + y))
   @test !equal(Ideal(R, x), Ideal(R, x, y))
   @test !isequal(Ideal(R, x, y, x - y), Ideal(R, x, y, x + y))
   @test isequal(Ideal(R, x, y), Ideal(R, x, y))
end

@testset "sideal.lead" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])
   i1 = @inferred lead(Ideal(R, x^2 + x*y + 1, 2y^2 + 3, R(7)))
   i2 = @inferred Ideal(R, x^2, 2y^2, R(7))
   @test isequal(i1, i2) || equal(i1, i2)

   r, (x, y) = polynomial_ring(Fp(7), ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x + y; 0 0]))
   i1 = @inferred lead(Ideal(R, x^2 + x*y + 1, 2y^2 + 3, R(7)))
   i2 = @inferred Ideal(R, x^2, 2y^2, R(7))
   @test isequal(i1, i2) || equal(i1, i2)

   R, (x, y) = FreeAlgebra(QQ, ["x", "y"], 10)
   i1 = @inferred lead(Ideal(R, x^2 + x*y + 1, 2y^2 + 3, R(7)))
   i2 = @inferred Ideal(R, x^2, 2y^2, R(7))
   @test isequal(i1, i2) || equal(i1, i2)
end

@testset "sideal.local" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:negdegrevlex)

   I = Ideal(R, y, x^2, (1 + y^3) * (x^2 - y))
   M = @inferred Singular.minimal_generating_set(I)

   @test equal(I, Ideal(R, x^2, y))
end

@testset "sideal.intersection" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I1 = @inferred Ideal(R, x^2 + x*y + 1, 2y^2 + 3, R(7))
   I2 = @inferred Ideal(R, x*y^2 + x + 1, 2x*y + 1, 7x + 1)
   I3 = @inferred Ideal(R, x, y)

   I = @inferred intersection(I1, I2)

   @test contains(I1, I)
   @test contains(I2, I)

   I = @inferred intersection(I1, I2, I3)

   @test equal(I, @inferred intersection(intersection(I1, I2), I3))
   @test equal(I, @inferred intersection(I1, intersection(I2, I3)))

   r, (x, y) = polynomial_ring(Fp(7), ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [1 1; 0 1]),
                           Singular.Matrix(r, [0 x + y; 0 0]))

   I1 = @inferred Ideal(R, x^2 + x*y + 1, 2y^2 + 3, R(7))
   I2 = @inferred Ideal(R, x*y^2 + x + 1, 2x*y + 1, 7x + 1)
   I3 = @inferred Ideal(R, x, y)

   I = @inferred intersection(I1, I2)

   @test contains(I1, I)
   @test contains(I2, I)
end

@testset "sideal.quotient" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = @inferred Ideal(R, x^2 + x*y + 1, 2y^2 + 3)
   J = @inferred Ideal(R, x*y^2 + x + 1, 2x*y + 1)
   K = @inferred Ideal(R, x*y + 1)

   A = @inferred quotient(I, J + K)
   B = @inferred intersection(quotient(I, J), quotient(I, K))

   @test equal(A, B)

   F, (q,) = FunctionField(QQ, ["q"])
   r, (x, y) = polynomial_ring(F, ["x", "y"], ordering = :degrevlex)
   R, (x, y) = GAlgebra(r, Singular.Matrix(r, [0 q; 0 0]),
                           Singular.Matrix(r, [0 0; 0 0]))

   I = std(Ideal(R, [x^3 + 2*x*y^2 + 2*x^2*y, y], twosided = true))
   J = std(Ideal(R, [x^2, x+ y], twosided = true))
   # the "equal" check here doesn't work because reduce only works with left ideals
   @test isequal(quotient(I, J), Ideal(R, [y, x^2]))
end

@testset "sideal.saturation" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Ideal(R, (x^2 + x*y + 1)*(2y^2+1)^3, (2y^2 + 3)*(2y^2+1)^2)
   J = Ideal(R, 2y^2 + 1)

   S, k = @inferred saturation(I, J)
   @test equal(S, Ideal(R, 2y^2 + 3, x^2 + x*y + 1))

   S2, k = @inferred saturation2(I, J)
   @test equal(S, S2)

   I = Ideal(R, (x*y + 1)*(2x^2*y^2 + x*y - 2) + 2x*y^2 + x, 2x*y + 1)
   J = Ideal(R, x)

   S, k = @inferred saturation(I, J)
   a = @inferred satstd(I, J)
   b = @inferred std(S)
   @test equal(a, b)

   S2, k = @inferred saturation2(I, J)
   @test equal(S, S2)

   r, (x, y) = polynomial_ring(ZZ, ["x", "y"], ordering = ordering_dp())
   I = Ideal(r, x^5 + x^4, y^2 + y)
   @test equal(satstd(I), Ideal(r, x + 1, y + 1))

end

@testset "sideal.slimgb" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Ideal(R, x^2 + x*y + 1, 2y^2 + 3)
   J = Ideal(R, 2*y^2 + 3, x^2 + x*y + 1)

   A = @inferred slimgb(I)

   leadA = @inferred lead(A)
   @test isequal(leadA, Ideal(R, 2y^2, x^2)) ||
         isequal(leadA, Ideal(R, x^2, 2*y^2))
   @test A.isGB == true

   B = slimgb(I, complete_reduction=true)

   @test isequal(B, Ideal(R, 2y^2 + 3, x^2 + x*y + 1)) ||
         isequal(B, Ideal(x^2 + x*y + 1, 2y^2 + 3))
   @test B.isGB == true
end

@testset "sideal.std" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Ideal(R, x^2 + x*y + 1, 2y^2 + 3)
   J = Ideal(R, 2*y^2 + 3, x^2 + x*y + 1)

   A = std(I)

   @test isequal(lead(A), Ideal(R, 2y^2, x^2)) ||
         isequal(lead(A), Ideal(R, x^2, 2*y^2))
   @test A.isGB == true

   B = std(I, complete_reduction=true)

   @test isequal(B, Ideal(R, 2y^2 + 3, x^2 + x*y + 1)) ||
         isequal(B, Ideal(x^2 + x*y + 1, 2y^2 + 3))
   @test B.isGB == true
end

@testset "sideal.interreduce" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

   I = @inferred Ideal(R, z*x + y^3, z + y^3, z + x*y)
   A = @inferred interreduce(I)
   @test equal(A, I)

   r, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering = :degrevlex)
   R, (x, y, z) = GAlgebra(r, Singular.Matrix(r, [0 1 1; 0 0 1; 0 0 0]),
                              Singular.Matrix(r, [0 x z; 0 0 y; 0 0 1]))

   I = @inferred Ideal(R, z*x + y^3, z + y^3, z + x*y)
   A = @inferred interreduce(I)
   @test equal(A, I)
end

@testset "sideal.fglm" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering = :lex)
   I = Ideal(R, y^3+x^2, x^2*y+x^2, x^3-x^2, z^4-x^2-y)
   J1 = @inferred std(I, complete_reduction = true)
   J2 = @inferred fglm(I, :degrevlex)

   @test gens(J1) == gens(J2)
end

@testset "sideal.reduction" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   f = x^2*y + 2y + 1
   g = y^2 + 1

   I = Ideal(R, (x^2 + 1)*f + (x + y)*g + x + 1, (2y^2 + x)*f + y)
   J = @inferred std(Ideal(R, f, g))

   @test isequal(reduce(I, J), Ideal(R, x + 1, y))

   h = (x^2 + 1)*f + (x + y)*g + x + 1

   @test (@inferred reduce(h, J)) == x + 1
end

@testset "sideal.free_resolution" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Ideal(R, x^2*y + 2y + 1, y^2 + 1)

   F1 = @inferred fres(std(I), 4)
   F2 = @inferred sres(std(I), 4)
   F3 = @inferred nres(I, 4)
   F4 = @inferred mres(I, 4)
   F5,TT = @inferred mres_with_map(I, 4)

   # check resolution is of the correct length
   @test (@inferred length(F1)) == 2
   @test (@inferred length(F2)) == 2
   @test (@inferred length(F3)) == 2
   @test (@inferred length(F4)) == 2

   # check index 1 is an ideal
   @test F1[1] isa sideal
   @test F2[1] isa sideal
   @test F3[1] isa sideal
   @test F4[1] isa sideal

   M1 = @inferred Singular.Matrix(F1[1])
   N1 = @inferred Singular.Matrix(F1[2])

   M2 = @inferred Singular.Matrix(F2[1])
   N2 = @inferred Singular.Matrix(F2[2])

   M3 = @inferred Singular.Matrix(F3[1])
   N3 = @inferred Singular.Matrix(F3[2])

   M4 = @inferred Singular.Matrix(F4[1])
   N4 = @inferred Singular.Matrix(F4[2])

   M5 = @inferred Singular.Matrix(F5[1])
   N5 = @inferred Singular.Matrix(F5[2])

   # check we have a complex
   @test iszero(M1*N1)
   @test iszero(M2*N2)
   @test iszero(M3*N3)
   @test iszero(M4*N4)
   @test iszero(M5*N5)
   @test iszero(Singular.Matrix(I)*TT-M5)
end

@testset "sideal.syzygy" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Ideal(R, x^2*y + 2y + 1, y^2 + 1)

   F = syz(I) #TODO @inferred broken

   M = @inferred Singular.Matrix(I)
   N = @inferred Singular.Matrix(F)

   # check they are actually syzygies
   @test iszero(M*N)
end

@testset "sideal.kernel" begin
   # twisted cubic
   P1, (t_0, t_1) = polynomial_ring(QQ, ["t_0", "t_1"])
   P3, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"])
   I = @inferred Ideal(P1, t_0^3, t_0^2*t_1, t_0*t_1^2, t_1^3)
   J = @inferred kernel(P3, I)
   @test ngens(J) == 3 && J[1] == z^2-y*w && J[2] == y*z-x*w && J[3] == y^2-x*z
end

@testset "sideal.eliminate" begin
   R, (x, y, t) = polynomial_ring(QQ, ["x", "y", "t"])

   I = Ideal(R, x - t^2, y - t^3)
   J = @inferred eliminate(I, t)
   @test equal(J, Ideal(R, x^3 - y^2))
   K = @inferred eliminate(I, t, x)
   @test equal(K, Ideal(R))

   r, (e, f, h, a) = polynomial_ring(QQ, ["e", "f", "h", "a"], ordering = :deglex)
   R, (e, f, h, a) = GAlgebra(r, Singular.Matrix(r, [0 1 1 1; 0 0 1 1; 0 0 0 1; 0 0 0 0]),
                                 Singular.Matrix(r, [0 -h 2*e 0; 0 0 -2*f 0; 0 0 0 0; 0 0 0 0]))

   I = @inferred Ideal(R, [e^3, f^3, h^3 - 4*h, 4*e*f + h^2 - 2*h - a])
   J = @inferred eliminate(I, e, f, h)
   @test equal(J, Ideal(R, [a^3 - 32*a^2 + 192*a]))
end

@testset "sideal.jet" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

   I = Ideal(R, x^5 - y^2, y^3 - x^6 + z^3)
   J = @inferred jet(I, 3)
   @test equal(J, Ideal(R, -y^2, y^3 + z^3))

   R, (x, y, dx, dy) = WeylAlgebra(QQ, ["x", "y"])

   I = Ideal(R, x^5 - y^2, y^3 - x^6 + dx^3)
   J = @inferred jet(I, 3)
   @test equal(J, Ideal(R, -y^2, y^3 + dx^3))
end

@testset "sideal.zerodim" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"]; ordering=:negdegrevlex)

   I = Ideal(R, 3*x^2 + y^3, x*y^2)

   I = @inferred std(I)

   dim = @inferred vdim(I)

   J = @inferred kbase(I)

   B = Ideal(R, R(1), x, y, x*y, y^2, y^3, y^4)

   f = @inferred highcorner(I)

   # Check dimension
   @test dim == 7

   # Check vector space basis
   @test equal(J, B)

   #Check highcorner
   @test f == y^4

   R, (x, y, z) = polynomial_ring(Fp(32003), ["x", "y", "z"])
   @test ngens(kbase(std(Ideal(R, x^2, y^3, x*y*z)), 2)) == 5
end

@testset "sideal.independent_set" begin
   R, (x, y, u, v, w) = polynomial_ring(QQ, ["x", "y", "u", "v", "w"])

   I = Ideal(R, x*y*w, y*v*w, u*y*w, x*v)

   I = @inferred std(I)

   L1 = maximal_independent_set(I)

   L2 = maximal_independent_set(I, all = true)

   L3 = @inferred independent_sets(I)

   @test L1 == [x, y, u]
   @test L2 == [[x, y, u], [y, u, v], [x, u, w], [u, v, w]]
   @test L3 == [[x, y, u], [y, u, v], [x, u, w], [u, v, w], [y, w]]
   @test typeof(L1) == Vector{spoly{n_Q}}
   @test typeof(L2) == Vector{Vector{spoly{n_Q}}}
end

@testset "sideal.lift" begin
  R,(x,y) = polynomial_ring(QQ, ["x", "y"], ordering = ordering_ds())
  f = x^5 + x^2*y^2 + y^5
  G = Ideal(R, derivative(f, x), derivative(f, y))

  for I in [Ideal(R, f), Ideal(R, f^2)]
    (T, rest, U) = division(I, G)
    @test Matrix(I)*Matrix(U) == Matrix(G)*Matrix(T) + Matrix(rest)

    (T, rest, U) = lift(G, I, true, false, true)
    @test Matrix(I)*U == Matrix(G)*Matrix(T) + Matrix(rest)
  end
end

@testset "sideal.lift_std" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Ideal(R, x, y)
   G, T = @inferred Singular.lift_std(I)
   @test G[1] == y
   @test G[2] == x
   @test T[1,1] == 0
   @test T[2,2] == 0
   @test T[2,1] == 1
   @test T[1,2] == 1

   G, T = @inferred Singular.lift_std(I, complete_reduction = true)
   @test G[1] == y
   @test G[2] == x
   @test T[1,1] == 0
   @test T[2,2] == 0
   @test T[2,1] == 1
   @test T[1,2] == 1

   G, T, S = Singular.lift_std_syz(I)
   mG = Singular.Matrix(G)
   mI = Singular.Matrix(I)
   mS = Singular.Matrix(S)
   @test iszero(mI*mS)
   # @test mI*T == mG
   # should be true, but there are additional (invisible) zeros in G
   # causing a dimension mismatch
end

@testset "sideal.hilbert_series" begin
   R, (x, y, z) = polynomial_ring(Fp(32003), ["x", "y", "z"])
   p = hilbert_series(std(Ideal(R, x^2, y^2, z^2)))
   pop!(p)
   @test p == [1,0,-3,0,3,0,-1]

   R, (x, y, z) = polynomial_ring(Fp(32003), ["x", "y", "z"])
   p = hilbert_series(std(Ideal(R, x, y, z)), [1, 2, 3])
   pop!(p)
   @test p == [1,-1,-1,0,1,1,-1]

   R, (x, y, z) = polynomial_ring(Fp(32003), ["x", "y", "z"], ordering = ordering_dp(3))
   i = Ideal(R, [1*x^2*y+151*x*y*z^10+169*y^21,
                 1*x*z^14+6*x^2*y^4+3*z^24,
                 5*y^10*z^10*x+2*y^20*z^10+y^10*z^20+11*x^3])
   h = hilbert_series(std(i))

   R, (x, y, z) = polynomial_ring(Fp(32003), ["x", "y", "z"], ordering = ordering_lp(3))
   j = Ideal(R, [1*x^2*y+151*x*y*z^10+169*y^21,
                 1*x*z^14+6*x^2*y^4+3*z^24,
                 5*y^10*z^10*x+2*y^20*z^10+y^10*z^20+11*x^3])
   @test ngens(std_hilbert(j, h, complete_reduction = true)) ==
         ngens(std(j, complete_reduction = true))


   R, (x, y, z) = polynomial_ring(Fp(32003), ["x", "y", "z"], ordering = ordering_dp(3))
   i = Ideal(R, [1*x^2*y+151*x*y*z^10+169*y^21,
                 1*x*z^14+6*x^2*y^4+3*z^24,
                 5*y^10*z^10*x+2*y^20*z^10+y^10*z^20+11*x^3])
   w = [10,1,1]
   h = hilbert_series(std(i), w)
   R, (x, y, z) = polynomial_ring(Fp(32003), ["x", "y", "z"], ordering = ordering_lp(3))
   j = Ideal(R, [1*x^2*y+151*x*y*z^10+169*y^21,
                 1*x*z^14+6*x^2*y^4+3*z^24,
                 5*y^10*z^10*x+2*y^20*z^10+y^10*z^20+11*x^3])

   @test ngens(std_hilbert(j, h, w, complete_reduction = true)) ==
         ngens(std(j, complete_reduction = true))
   Qt,(t,)= polynomial_ring(QQ, ["t"])
   I=Ideal(R,[x,y,z])
   I=std(I)
   @test hilbert_series(I,Qt) == -t^3+3*t^2-3*t+1
   A,x = polynomial_ring(QQ,["x$i" for i in 1:37])
   I = Ideal(A,[2*x[11] - 2*x[17] - 2*x[24] + 2*x[32] - 111916*x[37], 2*x[4] - 2*x[8] - 2*x[26] + 2*x[34] - 41216*x[37], 2*x[2] - 2*x[9] - 2*x[20] + 2*x[35] + 37974*x[37], x[28] - x[36], x[21] - x[36], x[27] - x[28] + x[33] + x[36], x[26] - x[27] - x[33] + x[34], x[20] - x[21] + x[35] + x[36], x[15] - x[21] - x[28] + x[36], x[10] - x[36], x[25] - x[28] + x[31] + x[36], x[24] - x[25] - x[26] + x[27] - x[31] + x[32] + x[33] - x[34], -x[14] + x[15] + x[18] - x[21] + x[25] - x[28] + x[31] + x[36], x[13] - x[14] + x[18] - x[19] - 2*x[20] + 2*x[21] - x[26] + x[27] + x[33] - x[34] - 2*x[35] - 2*x[36], x[9] - x[10] + x[35] + x[36], x[6] - x[10] - x[28] + x[36], x[19] - x[21] + x[30] + x[36], -x[18] + x[19] + x[23] - x[25] - x[27] + x[28] + x[30] - x[31] - x[33] - x[36], x[17] - x[19] - x[30] + x[32], x[12] - x[14] - x[17] + x[18] - x[27] + x[28] + x[31] - x[32] - x[33] - x[36], x[8] - x[10] + x[34] + x[36], x[5] - x[6] - x[8] + x[10] - x[27] + x[28] - x[34] - x[36], x[3] - x[10] - x[21] + x[36], -x[18] + x[19] + x[20] - x[21] + x[29] + x[30] + x[35] + x[36], x[22] + x[23] + x[24] - x[25] - x[29] - x[30] - x[31] + x[32], x[16] + x[17] + x[18] - x[19] - x[22] - x[23] - x[24] + x[25], x[11] + x[12] + x[13] - x[14] - x[16] - x[17] - x[18] + x[19] + x[22] + x[23] + x[24] - x[25] + x[29] + x[30] + x[31] - x[32], x[7] + x[8] + x[9] - x[10] - x[33] + x[34] + x[35] + x[36], x[4] + x[5] + x[9] - x[10] + x[26] - x[27] + x[35] + x[36], x[2] + x[3] + x[9] - x[10] + x[20] - x[21] + x[35] + x[36], x[1] - x[3] - x[6] + x[10] - x[15] + x[21] + x[28] - x[36], -x[27]*x[36] + x[34]*x[35], -x[25]*x[36] + x[32]*x[35], x[14]*x[36] + x[19]*x[35] + x[25]*x[36] + x[27]*x[36] - x[32]*x[35] - x[34]*x[35], -x[19]*x[36] - x[25]*x[36] + x[32]*x[34] + x[32]*x[35], -x[19]*x[35] - x[19]*x[36] + x[25]*x[34] - x[25]*x[36] + x[32]*x[34] + x[32]*x[35], x[14]*x[36] - x[19]*x[35] + x[25]*x[34] + x[27]*x[32], x[14]*x[35] - x[14]*x[36] + x[19]*x[35] - x[19]*x[36] + x[25]*x[27] - x[25]*x[34] - x[27]*x[32] + x[32]*x[34], x[14]*x[34] + x[19]*x[27] - 2*x[19]*x[35] + 2*x[25]*x[34] - x[25]*x[36] + x[32]*x[35], x[14]*x[32] - 2*x[14]*x[36] + x[19]*x[25] - 2*x[19]*x[35] - x[27]*x[36] + x[34]*x[35]])
   I=std(I)
   @test hilbert_series(I,Qt)==-t^37 + 31*t^36 - 456*t^35 + 4200*t^34 - 26775*t^33 + 121737*t^32 - 376992*t^31 + 556512*t^30 + 1739100*t^29 - 16811300*t^28 + 75314624*t^27 - 246484224*t^26 + 650872404*t^25 - 1444243500*t^24 + 2750940000*t^23 - 4555556640*t^22 + 6611884290*t^21 - 8454204990*t^20 + 9552774000*t^19 - 9552774000*t^18 + 8454204990*t^17 - 6611884290*t^16 + 4555556640*t^15 - 2750940000*t^14 + 1444243500*t^13 - 650872404*t^12 + 246484224*t^11 - 75314624*t^10 + 16811300*t^9 - 1739100*t^8 - 556512*t^7 + 376992*t^6 - 121737*t^5 + 26775*t^4 - 4200*t^3 + 456*t^2 - 31*t + 1
   v=hilbert_series_data(I)
   @test v[1] == BigInt(-1)
   @test length(v) == total_degree(hilbert_series(I,Qt))+1
end

@testset "sideal.oscar#1702" begin
   F, (b2, b3, d1, d2, d3, m1, m2, e, k, k2, k3) = FunctionField(QQ, ["b2", "b3", "d1", "d2", "d3", "m1", "m2", "e", "k", "k2", "k3"]);
   R, (H, S, C) = polynomial_ring(F, ["H", "S", "C"]; ordering = ordering_lp());
   I = Ideal(R,-e*H*S - m1*k*H*C + m1*C, -e*H*S - m2*k2*S*C + m2*C, e*H*S - b3*k3*C^2 + b3*C);
   @test dimension(std(I)) == 1
end

@testset "sideal.divrem and reduce" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])
   a = Ideal(R,[x^3+y^3+x*y])
   b = Ideal(R, [x])
   q,r,u = divrem(a, b)
   @test Singular.Matrix(a)*Singular.Matrix(u) == Singular.Matrix(b)*Singular.Matrix(q)+Singular.Matrix(r)
   @test gens(reduce(a, b, complete_reduction=false))[1] ==y^3+x*y
   @test gens(reduce(a, b, complete_reduction=true))[1] == y^3
   q2,r2,u2 = divrem2(a, b)
   @test Singular.Matrix(a)*Singular.Matrix(u2) == Singular.Matrix(b)*Singular.Matrix(q2)+Singular.Matrix(r2)
   @test isequal(r2, reduce(a, b))
end
