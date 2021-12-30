include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/bechmarkingEveryProcedure/GroebnerWalkFinalBenchmarkProcedures.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/bechmarkingEveryProcedure/runbenchmark.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkingAlg/GroebnerWalkFinalBenchmark.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkingAlg/runbenchmark2.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkHelper")
include("readWriteHelper.jl")
using DataFrames
using CSV

function benchmarkVar3()
    cd("/Users/JordiWelp/Results2")
    prepare()
    prepare2()
    prepareAlloc()
    dim = 3
    ve = [1, 1, 1]
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (a, b, c) = Singular.PolynomialRing(
        Singular.QQ,
        ["a", "b", "c"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)
    ideals = []

    push!(
    ideals,
    Singular.Ideal(
        R,[48*b+23*a*c+47*a^5+69*a*b^4,
53*b*c^2,
65+21*a*b*c^2]))
push!(
    ideals,
    Singular.Ideal(
        R,[R(47),
4*a*b^3+64*a^3*b^2,
40*a+90*b^2+59*a*b^2+74*b^4*c]))
push!(
    ideals,
    Singular.Ideal(
        R,[89*c+78*a*c^3,
26*a*c^2,
26+4*c^2+89*a^4*b+73*a^2*b^2*c]))
for i in ideals
    runAll("sparseid(3,0,5,100-(1*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[65+a^3+37*a^4+60*a^3*c^2,
47*a*b+4*a^2*b+35*b*c^4,
a+23*a^3*b+84*a*b^2*c^2]))
push!(
    ideals,
    Singular.Ideal(
        R,[86+77*a+13*a^3+23*a^2*b+21*c^4+4*a*b^4+7*c^5,
82*c^2,
23*c^4+4*a^3*b*c]))
push!(
    ideals,
    Singular.Ideal(
        R,[71*c^3+67*a^4+14*b*c^4,
64+87*b^2+40*a^2*b+21*b^2*c^2+73*a^4*b,
50*b+21*a^2*b^3]))
for i in ideals
    runAll("sparseid(3,0,5,100-(2*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[22+38*b*c^2+60*a*c^4,
14*a+67*a^2+58*b^2*c^2+87*b^5+88*a^2*b*c^2,
67*a*c+83*a^2*c+87*b^4+18*a*c^3+74*b^3*c^2]))
push!(
    ideals,
    Singular.Ideal(
        R,[49*a+88*a*b*c^2+55*a*b*c^3,
5*c^2+83*a^3*b+10*a^4*c+28*a*c^4,
54+19*b^2+8*b*c^2+54*c^3+10*a*b^2*c+40*c^5]))
push!(
    ideals,
    Singular.Ideal(
        R,[44*c+60*b*c+3*a^2*c+43*c^4,
85*a^2+36*b^3*c^2+19*c^5,
5+64*a*b^2+55*a^2*b^2+18*b^4+26*b*c^4+5*c^5]))
for i in ideals
    runAll("sparseid(3,0,5,100-(3*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[74*a*c+48*b*c+26*a*b^2*c+9*c^4+83*a*b^4+45*b^5+49*a^3*c^2,
85+42*a*b^2*c+68*a*c^3+50*b^2*c^3+45*c^5,
26*a+59*a^2*b+49*a^2*c+75*c^3+40*a^2*c^3]))
push!(
    ideals,
    Singular.Ideal(
        R,[48*a*c+19*a*c^2+23*a^3*b^2+42*a*b^2*c^2,
18*b^2+32*a^2*b+64*b*c^2+44*a*b^3+2*b^4+23*a^3*b^2+2*a^2*b^3+13*b^3*c^2,
41+19*c+51*a*b^3+9*b*c^3+50*a*b*c^3]))
push!(
    ideals,
    Singular.Ideal(
        R,[49*c^3+78*c^5,
12*a+22*b^4+40*a^3*b*c+44*a*b*c^3,
75+87*a*b+86*c^2+16*a^2*b+53*c^3+87*a^3*c+87*a*c^3+28*b*c^3+23*a^2*b^2*c+15*b^4*c+65*b^2*c^3]))
for i in ideals
    runAll("sparseid(3,0,5,100-(4*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[36*b*c+80*b^3+29*c^3+70*a*b^2*c+23*a^2*c^2+86*a^3*b^2+29*b^4*c+55*b*c^4,
65+75*a*c+58*b^3+69*a^3*b+64*a*b^2*c+82*a*b^4+24*c^5,
53*b+39*a^2*b*c+42*a^3*b*c+42*a^2*b*c^2]))
push!(
    ideals,
    Singular.Ideal(
        R,[78+7*a*c+21*a^2*b+20*b^2*c+85*b^3*c+34*b*c^3+9*b^2*c^3,
37*a^2*b^2+89*a^3*b^2+45*a*b^4+76*a^3*c^2+21*a*c^4+22*b*c^4,
37*c+24*b^2+49*b^2*c+15*a^2*b*c+55*a*c^3+5*a^2*c^3]))
push!(
    ideals,
    Singular.Ideal(
        R,[18+89*b^2+80*a^3+22*a^2*c^2+75*b^2*c^2+58*b*c^3+64*a*b^3*c+85*a^3*c^2+4*a^2*b*c^2+46*a^2*c^3,
16*a*b+31*a^3+10*b^3+9*a*b^2*c+67*a^2*c^2+46*a^2*b*c^2,
86*a+19*a^5+8*b^4*c]))
for i in ideals
    runAll("sparseid(3,0,5,100-(5*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[48*a+18*a^2+90*a*c+7*b*c+81*a*b*c+2*b^2*c+69*b^4+57*a*b*c^2,
20+35*b+24*b^4+11*b^3*c+49*a*c^3+22*a^3*b^2+76*a*b^4+3*a^2*b^2*c+47*b^2*c^3,
89*a^2*b+9*a*b^2+62*a*b^3+74*a^4*b+32*a^3*b^2+87*a^2*b^3+62*b*c^4]))
push!(
    ideals,
    Singular.Ideal(
        R,[48+15*a+42*a^2+35*c^2+72*b^4,
37*c^2+17*a*b^2+46*a*b*c+78*a^2*b*c+33*a*b*c^2+34*b^2*c^2+63*a*c^3+71*a^3*b^2+22*a^4*c+29*a^2*b^2*c+12*a^3*c^2+66*a*c^4+86*b*c^4,
73*c+52*a*b^2+69*a*b*c+18*a*b*c^2+49*a^4*b+77*b^3*c^2]))
push!(
    ideals,
    Singular.Ideal(
        R,[50*b+63*c^2+51*a^3+26*b^4*c+77*a^3*c^2+67*b^3*c^2,
11*b+48*a*b+57*a^4+10*b^4+23*a^5+14*a*b^4,
19+79*a*c+3*a*b^2+27*a*b*c+74*b^2*c+62*a^3*c+31*a^2*b*c+30*a*b^2*c+43*b^3*c+8*a*b^3*c+31*a^2*b*c^2+86*b^3*c^2]))
for i in ideals
    runAll("sparseid(3,0,5,100-(6*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[42+23*a*b+8*a^2*b+20*a^2*c+76*a^3*b*c+78*b^3*c^2+2*a*b*c^3+6*b^2*c^3+25*b*c^4,
31*c+50*b*c+80*a^2*b+79*a^3*b+85*a*b*c^2+58*a*b^4+2*a^3*c^2+16*c^5,
88*a+52*a*b+17*b^3+78*b*c^2+40*a^2*b^2+14*b^3*c+85*a^2*c^2+85*b^2*c^2+80*c^4+19*a^2*b^2*c]))
push!(
    ideals,
    Singular.Ideal(
        R,[61*a*b+35*b*c^2+13*a*b*c^2+69*b^2*c^2+20*b*c^3+21*a*b^4+21*a*b^3*c+53*b^3*c^2,
90+39*a+62*a^3*c+61*b*c^3+55*c^4+22*b^4*c+35*a^3*c^2+41*a*c^4+75*c^5,
56*c+a*c+13*b*c+33*a*b^2+65*b^3+45*a^2*c+71*b*c^2+4*a*b*c^2+82*a*b^3*c+70*b^2*c^3]))
push!(
    ideals,
    Singular.Ideal(
        R,[87*a*c+73*a*b^2+62*a*c^2+55*a^2*b^2+78*b^4+77*c^4+2*a^3*b^2+40*b^2*c^3+49*b*c^4,
58+42*c+85*a*b+58*a^2*b*c,
70*c+11*b^2+20*a^2*b+47*b^2*c+69*b*c^2+4*a^3*b+72*a*b*c^2+72*a*c^3+43*a^2*b^3+76*a*b^4+22*a^3*b*c+36*a^2*b*c^2+64*b^2*c^3+28*c^5]))
for i in ideals
    runAll("sparseid(3,0,5,100-(7*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[29+47*c+89*b^2+76*b*c+26*a^3+5*b^2*c^2+66*c^4+21*a^4*b+31*a^3*b^2+33*a*b^3*c+18*a*b^2*c^2+51*b^2*c^3,
22*a+57*a*c^2+56*c^3+62*a^4+10*a^3*c+22*a*c^3+68*a^3*b*c+10*a^2*c^3+34*a*b*c^3,
21*b^2+15*a*b*c+42*b^2*c+83*a^3*b+81*a*b^2*c+80*b^2*c^2+33*a*b*c^3+42*a*c^4+82*b*c^4]))
push!(
    ideals,
    Singular.Ideal(
        R,[17+21*c+67*a^2*c+70*b^2*c+84*a^2*b^2+14*a*b^2*c+66*a*b*c^2+23*a*b^4+6*b^5+22*a^2*c^3,
b+56*a^2+7*b^2+42*a^2*c+36*a*b*c+28*c^3+53*a^4+49*a^3*c+12*a^2*b^3+53*a*b^4+6*b^4*c+30*a^2*c^3+88*a*b*c^3+54*b^2*c^3,
28*b*c+16*a^3*b+65*a^2*c^2+23*a*c^3+54*b^3*c^2+38*a*b*c^3]))
push!(
    ideals,
    Singular.Ideal(
        R,[64*b^2+4*a^2*b+78*c^3+17*a*b^3+48*b^4+28*a*b^4+16*a^3*c^2,
49*c+25*a^2+13*c^2+87*a^2*c+34*c^3+88*a^2*b*c+8*a^3*b^2+6*a^2*b*c^2+31*a*c^4,
64+30*b+14*a*b*c+47*a^2*b^2+38*a*b^3+23*a^3*c+10*a^2*b*c+21*b^3*c+46*a^4*c+38*a^3*b*c+34*b^4*c+17*a^3*c^2+89*a^2*b*c^2+71*b^3*c^2]))
for i in ideals
    runAll("sparseid(3,0,5,100-(8*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[88*b+63*c^2+49*a*b*c+78*b*c^2+31*a^3*b+77*a*b^2*c+12*b^3*c+44*a^2*c^2+49*b^2*c^2+36*a*c^3+42*a^3*b*c+21*a^2*b^2*c+44*a*b^2*c^2,
2+55*a*b+75*a*c+75*b*c+51*a*b*c+59*a^3*b*c+39*a^2*c^3+89*c^5,
36*c+35*a^2*b+68*b^2*c+2*c^3+55*a^3*c+20*a*b^2*c+77*a^2*c^2+30*a^5+61*a^4*b+42*a^4*c+29*a*b^3*c+77*a^3*c^2+69*c^5]))
push!(
    ideals,
    Singular.Ideal(
        R,[44+4*a+60*c^2+35*a^3+15*a^2*c+47*b^2*c+68*b*c^2+2*a*b^3+4*a^5+7*a^3*b*c+70*a^2*b^2*c+56*a*b^3*c+70*b^2*c^3+2*a*c^4+40*c^5,
31*c+34*b*c+52*a*c^2+33*a^3*b+45*a^3*c+82*a*b*c^2+77*b*c^3+30*a^2*b^2*c+31*a^2*b*c^2+83*c^5,
33*a*c+40*b*c+25*b^2*c+78*a^4+68*a^2*b^2+65*a*b*c^2+c^4+40*a^2*c^3+49*a*c^4]))
push!(
    ideals,
    Singular.Ideal(
        R,[21*c^2+42*a^2*b+12*a*b^2+15*a*c^2+43*b*c^2+34*a*b^3+10*a^3*c+57*a^2*b*c+36*b^2*c^2+64*b^5+86*a^3*b*c+7*b^4*c+78*a*b*c^3+52*b^2*c^3+82*c^5,
45*a*c+63*a^3+8*a^4+80*b^2*c^2+76*a^4*b+46*a*b^2*c^2+9*b^3*c^2+11*a*c^4,
54+73*a+24*c+70*a^2+7*a*b+15*a^2*c+89*a^4+57*b^4+26*b*c^3+38*a^2*b^2*c+81*b^4*c]))
for i in ideals
    runAll("sparseid(3,0,5,100-(9*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]
push!(
    ideals,
    Singular.Ideal(
        R,[8*c+50*b^2+29*b*c+55*a^3+8*a^2*c+27*b^2*c+67*a*c^2+55*b^3*c+26*a^4*b+85*a^2*b^3+10*b^5+15*a*b^2*c^2+64*a^2*c^3,
76*a+64*a*b+20*b^2+24*a^3*b+79*b^4+11*b^2*c^2+5*a*c^3+88*b*c^3+3*c^4+67*a^3*b^2+37*a^2*b^3+47*a*b^3*c+51*a^2*c^3,
46+35*a*b^2+26*a*c^2+27*a^4+54*a*c^3+31*a^5+42*a^2*b^3+74*b^2*c^3+33*b*c^4]))
push!(
    ideals,
    Singular.Ideal(
        R,[62+85*b+37*c+9*a*c+56*a*c^2+76*b*c^2+55*c^3+70*a^4+20*a*b^3+21*a*c^3+73*a*b^2*c^2+69*a*c^4,
36*b^3+30*b^2*c+61*a^3*b+54*a^2*b^2+10*a*b^3+87*b^5+49*a^3*c^2+49*a*b^2*c^2+76*a^2*c^3+11*a*b*c^3+49*b^2*c^3,
74*b^2+12*a*c+11*b*c+40*a^2*c+63*a^4+86*b^4+17*a^2*c^2+86*a^3*b^2+88*b^5+56*a^2*b^2*c+13*a*b^2*c^2+22*c^5]))
push!(
    ideals,
    Singular.Ideal(
        R,[30*a+86*a*b*c+74*b^2*c+9*c^3+57*a^4+65*a^3*b+46*a*b^3+23*b^3*c+13*a*c^3+65*b*c^3+3*b^4*c+15*a^3*c^2+69*b^3*c^2+86*b*c^4,
1+57*a^2+7*a*b+37*a^3+15*a^2*c+79*a*b^3+80*b*c^3+42*a^3*b^2+38*a^2*b*c^2+15*b^2*c^3+32*a*c^4,
77*b+76*b^2+2*c^2+15*a*c^2+9*a^4+40*a^3*b^2+8*a^2*b^2*c+16*a^3*c^2+76*a^2*c^3+24*a*b*c^3]))
for i in ideals
    runAll("sparseid(3,0,5,100-(10*2),90)", i, S, StartOrd, TarOrd)
end
ideals =[]

end
