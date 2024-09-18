@testset "smodule.constructors" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)

   M = Singular.Module(R, v1, v2)

   S = parent(M)

   @test elem_type(S) == smodule{spoly{n_Q}}
   @test elem_type(ModuleClass{spoly{n_Q}}) == smodule{spoly{n_Q}}
   @test parent_type(smodule{spoly{n_Q}}) == ModuleClass{spoly{n_Q}}

   @test base_ring(S) == R
   @test base_ring(M) == R

   @test S isa AbstractAlgebra.Set
   @test M isa AbstractAlgebra.Module

   @test isa(M, smodule)
end

@testset "smodule.jet" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^5 + 1, 2x^3 + 3y^2, x^2)

   w1 = vector(R, x+1, x*y+1, y)
   w2 = vector(R, R(1), 2*x^3 + 3*y^2, x^2)

   M = Singular.Module(R, v1, v2)

   N = jet(M,3)

   P = Singular.Module(R, w1, w2)

   @test P[1] == N[1]
   @test P[2] == N[2]
end

@testset "smodule.local" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:negdegrevlex)

   v1 = vector(R, x, y^2)
   v2 = vector(R, y - x, y - y^2)
   v3 = v1 + v2

   w1 = vector(R, y, y)
   w2 = vector(R, x, y^2)

   M = Singular.Module(R, v1, v2, v3)
   MM = Singular.minimal_generating_set(M)

   @test MM[1] == w1 && MM[2] == -v2
end

@testset "smodule.manipulation" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)

   M = Singular.Module(R, v1, v2)

   @test rank(M) == 3

   @test ngens(M) == 2

   @test M[1] == v1
   @test M[2] == v2

   N = deepcopy(M)

   @test ngens(N) == 2

   @test N[1] == v1
   @test N[2] == v2
end

@testset "smodule.lead" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering = ordering_c()*ordering_ds())

   v = vector(R, 2*x^10, 2*x^2 + 3*y + 4*z^3, R(0))
   @test lead(v) == vector(R, 2*x^10, R(0), R(0))

   m = Singular.Module(R, v, vector(R, R(0), R(0), 2+x))
   n = lead(m)
   @test n[1] == vector(R, 2*x^10, R(0), R(0))
   @test n[2] == vector(R, R(0), R(0), R(2))
end

#= @testset "smodule.slimgb" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)
   v3 = x*v1 + y*v2 + vector(R, x, y + 1, y^2)

   M = Singular.Module(R, v1, v2, v3)

   G = slimgb(M; complete_reduction=true)

   @test ngens(G) == 3

   @test G[1] == vector(R, x, y + 1, y^2)
   @test G[2] == vector(R, x + 1, x*y + 1, y)
   @test G[3] == vector(R, x^2 + 1, 2*x + 3*y, x)

   @test G.isGB == true

   # Simply test the interfacing works in this case
   G2 = slimgb(M)

   @test G2.isGB == true
end =#

@testset "smodule.std" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, x + 1, x*y + 1, y)
   v2 = vector(R, x^2 + 1, 2x + 3y, x)
   v3 = x*v1 + y*v2 + vector(R, x, y + 1, y^2)

   M = Singular.Module(R, v1, v2, v3)

   G = std(M; complete_reduction=true)

   @test ngens(G) == 3

   @test G[1] == vector(R, x, y + 1, y^2)
   @test G[2] == vector(R, x + 1, x*y + 1, y)
   @test G[3] == vector(R, x^2 + 1, 2*x + 3*y, x)

   @test G.isGB == true

   # Simply test the interfacing works in this case
   G2 = std(M)

   @test G2.isGB == true
end

@testset "smodule.syz" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, (x + 1)*y, (x*y + 1)*y, y)
   v2 = vector(R, (x + 1)*x, (x*y + 1)*x, x)

   M = Singular.Module(R, v1, v2)

   Z = syz(M)

   @test ngens(Z) == 1

   @test Z[1] == vector(R, x, -y)
end

@testset "smodule.modulo" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, x)
   v2 = vector(R, y)

   A = Singular.Module(R, v1, v2)
   B = Singular.Module(R, v1)

   M = modulo(A,B)

   @test ngens(M) == 2

   @test M[1] == vector(R, R(1), R(0))
   @test M[2] == vector(R, R(0), x)
end

@testset "smodule.lift" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, x)
   v2 = vector(R, y)

   A = Singular.Module(R, v1, v2)
   B = Singular.Module(R, v1)

   M,r = lift(A,B)

   @test M[1] == vector(R,R(1),R(0))
   @test iszero(r[1])
end

@testset "smodule.eliminate" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   v1 = vector(R, x)
   v2 = vector(R, y)

   A = Singular.Module(R, v1, v2)

   M = eliminate(A,x)

   @test M[1] == v2
end

@testset "smodule.sres" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

   I = Singular.Ideal(R, y*z + z^2, y^2 + x*z, x*y + z^2, z^3, x*z^2, x^2*z)
   M = std(syz(I))

   F = sres(M, 0)
   F2 = mres(M,0)
   F3 = nres(M,0)

   # We have R^6 <- R^9 <- R^5 <- R^1
   # All references agree that when written as follows:
   # 0 <- M <- R^6 <- R^9 <- R^5 <- R^1 <- 0
   # the length is 3, as numbering starts at index 0 with R^6

   @test length(F) == 3
   @test F[1] isa smodule
   @test length(F2) == 2
   @test F2[1] isa smodule
   @test length(F3) == 2
   @test F3[1] isa smodule

   M1 = Singular.Matrix(F[1])
   M2 = Singular.Matrix(F[2])
   M3 = Singular.Matrix(F[3])

   @test Singular.Matrix(Singular.Module(M1)) == M1

   @test iszero(M1*M2)
   @test iszero(M2*M3)

   M1 = Singular.Matrix(F2[1])
   M2 = Singular.Matrix(F2[2])

   @test iszero(M1*M2)

   M1 = Singular.Matrix(F3[1])
   M2 = Singular.Matrix(F3[2])

   @test iszero(M1*M2)

   F = sres(M, 1)

   @test length(F) >= 1 # Singular can return more than asked for

   M1 = Singular.Matrix(F[1])
   M2 = Singular.Matrix(F[2])

   @test iszero(M1*M2)

   F = sres(M, 2)

   @test length(F) >= 2 # Singular can return more than asked for

   M1 = Singular.Matrix(F[1])
   M2 = Singular.Matrix(F[2])
   M3 = Singular.Matrix(F[3])

   @test iszero(M1*M2)
   @test iszero(M2*M3)

   F = sres(M, 4)

   @test length(F) == 3
end


@testset "smodule.vdim" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
   v1 = vector(R, x, R(0))
   v2 = vector(R, R(0), x)
   v3 = vector(R, y, R(0))
   v4 = vector(R, R(0), y)
   @test vdim(std(Singular.Module(R, v1, v2, v3, v4))) == -1

   R, (x, y) = polynomial_ring(QQ, ["x", "y"])
   v1 = vector(R, x, R(0))
   v2 = vector(R, R(0), x)
   v3 = vector(R, y, R(0))
   v4 = vector(R, R(0), y)
   @test vdim(std(Singular.Module(R, v1, v2, v3, v4))) == 2
end

@testset "smodule.hilbert_series" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
   v1 = vector(R, x^2, R(0))
   v2 = vector(R, y^2, R(0))
   v3 = vector(R, z^2, R(0))
   v4 = vector(R, R(0), R(1))
   M = Singular.Module(R, v1, v2, v3, v4)
   M.isGB = true
   @test hilbert_series(M,[1,1,1],[0,0]) == [1,0,-3,0,3,0,-1,0]
end

@testset "smodule.divrem and reduce" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])
   a = Singular.Module(R,vector(R, x^3+y^3+x*y, R(1)))
   b = Singular.Module(R, vector(R, x, R(0)))
   q,r,u = divrem(a, b)
   @test Singular.Matrix(a)*Singular.Matrix(u) == Singular.Matrix(b)*Singular.Matrix(q)+Singular.Matrix(r)
   r = reduce(a, b, complete_reduction=false)
   @test r[1] == vector(R, y^3+x*y, R(1))
   r = reduce(a, b, complete_reduction=true)
   @test r[1] == vector(R, y^3, R(1))
end

@testset "smodule.prune" begin
   R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

   v1 = vector(R, R(1),R(0),R(0))
   v2 = vector(R, x,R(0),R(0))
   v3 = vector(R, R(0),z,R(0))
   v4 = vector(R, R(0),R(0),R(1))

   M = Singular.Module(R, v1, v2, v3, v4)

   G,T = prune_with_map(M)

   @test ngens(G) == 1

   @test G[1] == vector(R, z)
   @test T[3,1] == R(1)

   G,T,p = prune_with_map_projection(M)
   @test p[1] == 1
   @test p[2] == 1
   @test p[3] == 2

end

@testset "smodule.contains" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Singular.Module(R, vector(R,y^3))
   J = Singular.Module(R, vector(R,y^2))

   @test contains(I, J) == false
   @test contains(J, I) == true
end

@testset "smodule.quotient" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Singular.Module(R, vector(R,x^2 + x*y + 1), vector(R,2y^2 + 3))
   J = Singular.Module(R, vector(R,x*y^2 + x + 1), vector(R,2x*y + 1), vector(R,x*y + 1))

   A = quotient(I, J)
   B = Singular.Module(R, vector(R, 2*y^2+3), vector(R, x^2+x*y+1))
   @test A[1] == B[1]
   @test A[2] == B[2]
end

@testset "smodule.saturation" begin
   R, (x, y) = polynomial_ring(QQ, ["x", "y"])

   I = Singular.Module(R, vector(R,(x^2 + x*y + 1)*(2y^2+1)^3), vector(R,(2y^2 + 3)*(2y^2+1)^2))
   J = Singular.Module(R, vector(R,2y^2 + 1))

   A,k = saturation(I, J)

   B = Singular.Module(R, vector(R,2*y^2+3), vector(R,x^2+x*y+1))
   @test A[1] == B[1]
   @test A[2] == B[2]
   @test k == 2

   A,k = saturation2(I, J)

   @test A[1] == B[1]
   @test A[2] == B[2]
   @test k == 2
end
