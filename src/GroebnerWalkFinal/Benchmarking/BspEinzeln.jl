include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/bechmarkingEveryProcedure/GroebnerWalkFinalBenchmarkProcedures.jl")
include("readWriteHelper.jl")
include("parser.jl")
include("runbenchmark.jl")


using DataFrames
using CSV
function runAllSingleExample()
    cd("/Users/JordiWelp/Results")
#Katsura 5
    dim = 5
    ve = [1, 1, 1, 1, 1]
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (v,u,t,z,y) = Singular.PolynomialRing(
        Singular.QQ,
        ["v", "u", "t", "z", "y"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 = 2 * y^2 + 2 * z^2 + 2 * t^2 + 2 * u^2 + v^2 - v
    f2 =  2*z*y + 2 * z * t + 2 * t * u + 2 * u * v - u
    f3 =  2 * y * t + 2 * z * u + 2*u^2 + 2 * t * v - t
    f4 =  2 * y * u + 2 * t * u + 2 * z * v - z
    f5 =  2 * y + 2 * z + 2 * t + 2 * u + v - 1
    I = Singular.std(Singular.Ideal(R, [f1, f2,f3,f4,f5]), complete_reduction=true)
    runb("Katsura5", I, S, StartOrd, TarOrd)

    s = Singular.std(
        Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(I)]),
        complete_reduction = true,
    )
    df = DataFrame(a = ["test1"], b = ["test2"],c=["example"])
    savew(df, "correct.txt")
    for id in ideals
        a = isequal(
            Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
            s,
        )
        b = equalitytest(s, id)
        df = DataFrame(a = [a], b = [b], c=["Katsura5"])
        savea(df, "correct.txt")
    end

    #Katsura6

    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (v,u,t,z,y,x) = Singular.PolynomialRing(
        Singular.QQ,
        ["v", "u", "t", "z", "y", "x"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)
    f1 = 2 * x^2 + 2 * y^2 + 2 * z^2 + 2 * t^2 + 2 * u^2 + v^2 - v
    f2 = 2*x * y +2* y * z + 2 * z * t + 2 * t * u + 2 * u * v - u
    f3 = 2 * x * z + 2 * y * t + 2 * z * u +2* u^2 + 2 * t * v - t
    f4 = 2 * x * t + 2 * y * u + 2 * t * u + 2 * z * v - z
    f5 = t^2 + 2 * x * v + 2 * y * v + 2 * z * v - y
    f6 = 2 * x + 2 * y + 2 * z + 2 * t + 2 * u + v - 1

    I = Singular.std(Singular.Ideal(R, [f1, f2,f3,f4,f5,f6]), complete_reduction= true)
    runb("Katsura6", I, S)

    s = Singular.std(
        Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(I)]),
        complete_reduction = true,
    )
    df = DataFrame(a = ["test1"], b = ["test2"],c=["example"])
    savew(df, "correct.txt")
    for id in ideals
        a = isequal(
            Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
            s,
        )
        b = equalitytest(s, id)
        df = DataFrame(a = [a], b = [b], c=["Katsura6"])
        savea(df, "correct.txt")
    end


        dim = 7
        ve = [1, 1, 1, 1, 1, 1, 1]
        example ="Cyclic7"
        StartOrd = ordering_as_matrix(ve, :lex)
        TarOrd = ordering_as_matrix(:lex, dim)
        R, (x1,x2,x3,x4,x5,x6,x7) = Singular.PolynomialRing(
            Singular.QQ,
            ["x1","x2","x3","x4","x5","x6","x7"],
            ordering = Singular.ordering_M(StartOrd),
        )
        S = change_order(R, TarOrd)
  f1=x1+x2+x3+x4+x5+x6+x7
  f2=x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x1*x7+x6*x7
  f3=x1*x2*x3+x2*x3*x4+x3*x4*x5+x4*x5*x6+x1*x2*x7+x1*x6*x7+x5*x6*x7
  f4=x1*x2*x3*x4+x2*x3*x4*x5+x3*x4*x5*x6+x1*x2*x3*x7+x1*x2*x6*x7+x1*x5*x6*x7+x4*x5*x6*x7
  f5=x1*x2*x3*x4*x5+x2*x3*x4*x5*x6+x1*x2*x3*x4*x7+x1*x2*x3*x6*x7+x1*x2*x5*x6*x7+x1*x4*x5*x6*x7+x3*x4*x5*x6*x7
  f6=x1*x2*x3*x4*x5*x6+x1*x2*x3*x4*x5*x7+x1*x2*x3*x4*x6*x7+x1*x2*x3*x5*x6*x7+x1*x2*x4*x5*x6*x7+x1*x3*x4*x5*x6*x7+x2*x3*x4*x5*x6*x7
  f7=x1*x2*x3*x4*x5*x6*x7-1

      I = Singular.std(Singular.Ideal(R, [f1, f2,f3,f4,f5,f6,f7]), complete_reduction= true)
      runb(example, I, S)

      s = Singular.std(
          Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(I)]),
          complete_reduction = true,
      )
      df = DataFrame(a = ["test1"], b = ["test2"],c=["example"])
      savew(df, "correct.txt")
      for id in ideals
          a = isequal(
              Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
              s,
          )
          b = equalitytest(s, id)
          df = DataFrame(a = [a], b = [b], c=[example])
          savea(df, "correct.txt")
end

dim = 6
ve = [1, 1, 1, 1, 1, 1]
example ="Cyclic6"
StartOrd = ordering_as_matrix(ve, :lex)
TarOrd = ordering_as_matrix(:lex, dim)
R, (x1,x2,x3,x4,x5,x6) = Singular.PolynomialRing(
    Singular.QQ,
    ["x1","x2","x3","x4","x5","x6"],
    ordering = Singular.ordering_M(StartOrd),
)
S = change_order(R, TarOrd)


  f1 =x1+x2+x3+x4+x5+x6
  f2 =x1*x2+x2*x3+x3*x4+x4*x5+x1*x6+x5*x6
  f3 =x1*x2*x3+x2*x3*x4+x3*x4*x5+x1*x2*x6+x1*x5*x6+x4*x5*x6
  f4 =x1*x2*x3*x4+x2*x3*x4*x5+x1*x2*x3*x6+x1*x2*x5*x6+x1*x4*x5*x6+x3*x4*x5*x6
  f5 =x1*x2*x3*x4*x5+x1*x2*x3*x4*x6+x1*x2*x3*x5*x6+x1*x2*x4*x5*x6+x1*x3*x4*x5*x6+x2*x3*x4*x5*x6
  f6 =x1*x2*x3*x4*x5*x6-1

  I = Singular.std(Singular.Ideal(R, [f1, f2,f3,f4,f5,f6]), complete_reduction= true)
  runb(example, I, S)

  s = Singular.std(
      Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(I)]),
      complete_reduction = true,
  )
  df = DataFrame(a = ["test1"], b = ["test2"],c=["example"])
  savew(df, "correct.txt")
  for id in ideals
      a = isequal(
          Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
          s,
      )
      b = equalitytest(s, id)
      df = DataFrame(a = [a], b = [b], c=[example])
      savea(df, "correct.txt")
end
end
