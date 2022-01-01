include("GroebnerWalkFinalBenchmarkProcedures.jl")
include("runbenchmark.jl")
include("GroebnerWalkFinalBenchmark.jl")
include("runbenchmark2.jl")
include("BenchmarkHelper")
include("readWriteHelper.jl")
include("Examples")

using DataFrames
using CSV

function benchmarkVar5()
    cd("/Users/JordiWelp/Results")
    prepare()
    prepare2()
    prepareAlloc()
    dim = 5
    ve = [1, 1, 1, 1, 1]
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (a, b, c, d, e) = Singular.PolynomialRing(
        Singular.QQ,
        ["a", "b", "c", "d", "e"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)
    ideals = []

    push!(
        ideals,
        Singular.Ideal(
            R,[5*b*e^4,
14+16*b+37*a*b*d+88*a*e^2+20*c^2*e^2+83*a^2*b*c*e+19*b^2*d^2*e+41*b*e^4,
50*b^2+41*b*c^3+20*a*c*e^2]))
    push!(
        ideals,
        Singular.Ideal(
            R,[49*c^3*d+46*a^3*b*d+78*d^5,
12+73*b*d^2+82*a*d*e+90*b*c*e^2+67*b^3*c^2+73*b^2*c^2*e,
20*a+87*b*c+55*a^3*c]))
    push!(
        ideals,
        Singular.Ideal(
            R,[23*d^2*e,
68+5*b*c+46*b*c*d+61*b^4+69*b^3*e+23*b^2*c^3+9*d^4*e,
64*a+74*a*c*d^2+4*a*c^3*d+4*c*d*e^3]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(1),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[45*c*d^2+78*a*b*c*d+16*c*d^2*e+24*b^5+55*a^3*b*d+57*b^2*d^2*e,
89*b+80*a*d+59*b*d^2+31*e^3+40*c^3*d+29*a*b*d*e^2,
35+47*a*b*d^2+55*c*d^2*e+13*b^3*c^2+27*a*b*c*d^2+84*b^2*c*e^2+43*b*c*d*e^2]))
push!(
    ideals,
    Singular.Ideal(
        R,[a^2*c+87*a^2*c^2+43*b^3*e+27*a*d^2*e+15*a^3*b*d+35*b*c*d^3+60*a*e^4,
67+40*b+11*a*d^2+59*a*b*e+50*b*e^4,
44*a*e+76*a^2*b^2+53*a*b^2*e+9*a*d^4+29*a*b*c*d*e+87*a*b*c*e^2+35*a*d*e^3]))
    push!(
        ideals,
        Singular.Ideal(
            R,[27+61*c^3*e+23*b*c*d*e+83*a*d^2*e+84*b^4*c+71*b*c^2*d^2+64*b^3*c*e+81*a*b*c^2*e,
42*b*e+82*a*b*d+47*a*b^2*d+39*b^2*e^2,
37*d+31*a^2*c+5*a*e^2+81*a^3*b*e+84*a^2*c*e^2+50*a*b*d*e^2+23*e^5]))

    for i in ideals
        runAll("sparseid(3,0,5,100-(2),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[55+90*c*d+64*d^2*e+16*d*e^2+77*a^2*b*c+90*b^3*d+73*a*b^2*e+c^5+69*b^2*c*d^2+42*a*c^2*d^2,
86*b^2+88*c*e^2+45*a^3*b+69*c^3*d+35*a*c^2*e+4*b^3*d^2+52*c*d^4+50*a^2*c^2*e+79*b^2*c*e^2+74*b*c^2*e^2,
78*d+61*a*c*d+23*b*d^3+7*a^3*c^2+70*b^2*c*e^2+77*a^2*d*e^2+78*c*d^2*e^2]))
    push!(
        ideals,
        Singular.Ideal(
            R,[44*c^2*d+15*b*c^3+80*a*d*e^2+87*a*c^3*d+14*a^4*e+51*a^3*c*e,
20*a*c*d+13*a^3*b+86*b*c^2*d+27*d^3*e+42*a*b*c^2*d+35*a*b^2*d*e+61*a*b*d*e^2+72*d^3*e^2+44*a*d*e^3+59*d*e^4,
84+9*d+50*b^2+2*b*d+44*b*c*d+77*e^3+46*d^3*e+a^2*e^2+87*b^3*c*d+57*a*b*c*e^2+78*a^2*e^3]))
    push!(
        ideals,
        Singular.Ideal(
            R,[83+19*d+7*a^2*c+26*b*d^2+61*d^3+2*a*d*e+83*a^2*c^2+49*b*d^3+3*a*c*d*e+33*a*b^2*d^2+55*a*b*d^2*e+35*a^2*c*e^2+45*a*d^2*e^2,
36*a^2*b*d+72*b*c^2*d+69*a*b^3*c+69*c^3*d^2+33*a^2*b^2*e+7*b^2*c^2*e+85*b*e^4,
37*b*e+80*e^2+66*a^2*c*e+46*b*d*e^2+79*c^3*d*e+31*b^2*d^2*e+82*b^2*e^3]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(3),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[69*b*d+51*a*b*c+46*a*c*d^2+28*b*c*d^2+34*b^3*e+73*b*c^2*e+74*d*e^3+61*e^4+79*a*b^2*c^2+3*b*c^2*d*e+16*a*c*d*e^2+52*b^2*e^3,
79+81*a+7*a*c^2+17*a^2*d+78*b^2*d+41*c*d^2+57*c^4+90*a^3*b^2+7*b*c*d^3+38*b^2*c^2*e+44*c^4*e,
53*d^2+81*b*c*d^2+16*a^2*c*e+48*a^4*c+55*a^2*b^2*c+18*a*b^2*c*d+9*d^5+87*b*c^2*d*e+62*a*d^3*e+7*a^2*d*e^2+87*a*b*e^3]))
    push!(
        ideals,
        Singular.Ideal(
            R,[90*c+56*a*e+81*c^2*e+36*b*d*e+55*a^2*b*d+60*a*b*e^2+37*a^2*b*c^2+65*b^4*d+36*a*b*c*e^2+46*d*e^4,
72+25*a*d+72*a*b^2*c+30*a^2*c^2+11*a^3*d+68*c^2*e^2+4*a*b^2*c*d+86*a^3*c*e+40*b^3*d*e+19*a*c*d^2*e+74*a*b*c*e^2,
87*a^2*b+27*b^2*c+41*d^3+86*b^4+3*a*c^2*e+75*a*d^2*e+22*a^4*c+11*a^3*b*c+89*a^2*c*d^2+18*b*c^2*d^2+61*a*c*d^3+4*b^3*d*e+38*b^2*d*e^2]))
    push!(
        ideals,
        Singular.Ideal(
            R,[36+11*a*b*d+6*b^2*e+74*a*b^2*c+85*a^3*e+24*a^2*e^2+40*c*d*e^2+42*a^4*c+5*a*b^2*c^2+67*a^3*c*d+33*a^2*c*d^2,
11*b*d+34*a*e+34*c*d^2+33*a*c*e+86*b*c*d*e+73*a^3*c^2+45*a*c*d^3+55*b^2*c^2*e+11*d^4*e+9*a^2*c*e^2,
56*c+30*b*c*e+48*a^2*c^2+82*c^3*d+75*b^2*d^2+23*c*e^3+51*a^2*b^3+71*a^4*d+14*a^2*c*d^2+29*c^2*d^3+54*a^3*b*e+36*a^3*d*e+4*a^2*d^2*e]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(4),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[46*d+22*a*b*c+71*a^2*e+18*a^4+69*b^3*d+47*a^2*c^3+37*c^5+88*a^3*d^2+13*a^2*c^2*e+88*a*c^2*d*e+19*b^2*c*e^2+48*a*b*d*e^2,
9+30*b*d+38*a^2*e+87*b^2*c*d+25*a*d^3+82*d^4+40*a^2*d^3+74*b*c^2*e^2+70*a*c*d*e^2+45*a*b*e^3,
74*a*c+46*e^2+19*a^2*b+12*b^3+64*a*c*d+76*b*c^3+90*b^2*d^2+40*c*d^3+38*b^3*e+35*b*c*e^2+5*a*e^3+79*b^4*c+44*a*b^3*e+69*a*c^3*e+39*a*b*d^2*e+35*b*c*d*e^2+9*b^2*e^3+71*b*d*e^3+54*e^5]))
    push!(
        ideals,
        Singular.Ideal(
            R,[84*e+34*a*b*c*d+69*b*c*d^2+70*a^2*e^2+9*a^5+47*a^2*b*c*d+63*a^2*d*e^2+45*b*c*d*e^2,
83+29*a*b+46*b^3+18*c*d*e+82*d^2*e+a^2*b^2+19*b^4+26*a*d^2*e+31*a*e^3+85*d*e^3+78*a^4*c+89*a*b^3*e+76*a*b*d*e^2+41*a*c*d*e^2+16*b*c*d*e^2+90*a*d^2*e^2+90*c*d^2*e^2+37*c^2*e^3+33*c*e^4,
85*a*c+83*c*e+17*a^3+82*a*b*d+65*c*e^2+82*a^3*b+86*a*b*c^2+77*b^2*d^2+19*a*b^4+27*b*c^4+47*a*b^2*c*d+86*b*c^2*d^2+25*b^3*c*e+53*a*c*e^3]))
    push!(
        ideals,
        Singular.Ideal(
            R,[13*b+44*a^2+39*a^2*b^2+12*a^2*b*e+56*d*e^3+41*e^4+3*a*b*c^3+85*b^4*d+42*a^2*d*e^2,
21*a*b+47*b*c^2+31*b*d^2+39*c*d^2+86*e^3+53*a*c^2*d+74*a*d^3+15*a*d*e^2+51*a^2*b^3+60*a^4*c+22*a*b*c^3+14*a*b*c*d^2+75*b*c*d^3+61*a*b^2*c*e+7*b^3*c*e+24*a*c*d^2*e+87*b*c*d*e^2+78*c^2*d*e^2+82*c*e^4,
13+68*d*e+26*a^2*d+23*b^2*d+12*a^3*d+40*a*b^2*e+5*b*c^2*e+57*b*e^3+60*a^2*b^2*d+b^3*d^2+31*a*c^3*e+61*b^2*d*e^2+4*a*b*e^3]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(5),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[70*b^2*d+45*a*b^2*c+90*b^3*c+49*a^2*c^2+22*c*d^3+81*b*d^2*e+62*b^5+68*a^3*b*c+86*a^2*b^2*c+70*a*b^3*c+15*a*c^4+21*a*b^3*d+20*a*b^2*c*d+44*b*c^3*d+33*a*b*d^3+17*a^2*b*c*e+53*d^4*e+90*c^2*e^3,
65*b+16*c*d+31*d*e+68*e^2+19*a^2*b+47*c^4+18*a*b*c*d+64*a^2*c*e+67*a*d^2*e+63*a^2*e^2+79*a^4*c+11*a^2*b^2*e+31*c^2*d^2*e+9*c*e^4,
31+26*a*b^2+86*a^2*d+71*a*b*e+26*a*c*e+63*b*e^2+87*a*c*d^2+59*b*c*d^2+23*a*e^3+18*a*b^2*c*d+88*b*c^2*d^2+33*c^2*d^3+28*a^4*e+5*b*c^2*e^2+17*b*d^2*e^2+86*c*e^4]))
    push!(
        ideals,
        Singular.Ideal(
            R,[90*a+56*a*c+3*a^2*c+53*b*d^2+90*b^4+60*a*b*c^2+38*c*d^3+66*a^2*d*e+4*b^2*d*e+11*a*e^3+34*b^3*c^2+16*b^3*c*d+59*a*b*c^2*d+12*a*b^2*d^2+23*b*d^4+29*a^3*d*e+90*c^2*d^2*e,
48+57*a*b+25*b^3+17*a^2*d+77*a*d^2+41*a*b*c*d+45*b*d^3+64*d^3*e+4*a*d*e^2+6*a*b^3*c+51*a*b^2*c^2+69*b^2*c^2*d+54*a*c^2*d^2+26*c^3*d^2+51*c*d^4+65*b*d^3*e+46*c*d^3*e+55*b^3*e^2+19*b*c^2*e^2+47*d*e^4,
30*b*d+40*a^2*b+69*c*d*e+15*a*b*c^2+23*c^3*d+85*d*e^3+64*a^3*c^2+75*b^3*c^2+39*c^5+49*a*c^3*d+49*c*d^4]))
    push!(
        ideals,
        Singular.Ideal(
            R,[63*b*e+76*e^2+60*a*c*d+76*e^3+22*a^3*c+72*a*c*e^2+83*c^2*e^2+33*a^2*b^2*c+45*a^2*b*c*d+90*b*c^3*e+3*a^3*d*e+15*c^3*d*e+75*a*c*e^3+28*b*c*e^3+77*c^2*e^3,
35+83*d*e+65*b*c^2+29*c*d^2+33*d^2*e+20*a^2*c^2+67*a^3*d+20*a^2*b*d+4*d^4+54*a*c*d*e+50*b*c*d*e+77*a^4*c+22*b^4*c+62*a*b*c^3+11*a^2*b^2*d+19*a^2*b*c*d+19*c^4*d+24*a*b*c*d^2+32*a*b*d^3+45*c*d^3*e+29*b^2*e^3+23*b*d*e^3,
45*a+74*b^2*c+17*d*e^2+26*a*b^2*d+22*c^3*d+2*c*e^3+17*d*e^3+46*a*b^3*c+73*b^3*c*d+3*a^2*b^2*e+4*c^4*e]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(6),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[4+25*d+56*b^2+48*b*d+22*a*b^2+88*b^2*d+37*a*b*e+80*c^4+10*b*c*d^2+60*c*d^2*e+39*c^2*e^2+60*a*b^3*c+53*a*c^4+29*c^5+32*a*b*c^2*d+8*a*c*d^3+21*a^2*b^2*e+14*a*b^3*e+68*c*e^4,
5*b*c+84*d^2+55*a*b^2+31*b^2*c+7*b^3*d+14*a^2*d^2+81*b^3*e+67*b*d^2*e+64*a*c*e^2+24*a^2*b^2*c+83*b^3*c^2+40*a^4*d+80*b^4*d+23*b*c^3*d+76*a^3*d^2+50*b*d^4+3*a^4*e+62*c^3*e^2+53*b^2*d*e^2+60*c^2*d*e^2+35*d^3*e^2+77*d*e^4,
86*b+81*c^2*d+72*b*d^2+63*d^3+10*a*b^3+14*b^2*c^2+64*b^3*d+47*a^3*e+61*a^2*c*e+25*c*d*e^2+48*c^4*d+88*a^2*c*d^2+58*a*c^2*d^2+20*b*c*d^3+82*b^3*c*e+85*b*d*e^3]))
    push!(
        ideals,
        Singular.Ideal(
            R,[7*c*d+81*a^3+48*a*c*e+32*a*b^2*e+75*a^2*c*e+59*a*c^2*e+41*a*b*d*e+57*d^2*e^2+8*b*e^3+33*c^4*d+45*a^2*b*d^2+40*b^4*e+39*a*d*e^3,
63+80*a+28*b+24*c*d+24*a^3*c+56*a*c^3+23*c^3*d+74*a^2*d^2+62*c*d^2*e+78*a^2*b^2*d+18*a^2*c*d^2+2*c*d^4+42*a^2*d^2*e+24*c*d^3*e+18*b^3*e^2+27*a*c^2*e^2+34*c^3*e^2+57*a^2*d*e^2+90*a*c*d*e^2+42*a*d*e^3+13*c*d*e^3,
75*a^2+60*b^2+68*c*d^2+88*d^3+86*a^2*e+4*a*b*e+65*c*d*e+73*e^3+71*b^3*d+65*b*c*d^2+85*c^3*e+21*c*d*e^2+29*a^2*b^3+20*a^2*b^2*c+90*a^2*b^2*d+5*a*b*d^3+26*a^2*b*c*e+87*a^3*d*e+25*a*b^2*d*e+90*a^2*c*e^2+34*b^2*e^3+55*c*d*e^3+83*a*e^4]))
    push!(
        ideals,
        Singular.Ideal(
            R,[25*a^2+13*c*d+33*d^3+28*a^2*e+22*a*b*e+29*a*b^3+62*a^2*c^2+37*a^3*d+5*b^2*c*d+64*b^2*c*e+59*a^2*d*e+54*a*d^2*e+62*d^3*e+8*a*b^3*c+43*a^2*c^3+64*a*b^2*c*e+28*a^2*c^2*e+39*c*d^2*e^2+26*d^3*e^2,
14*a*c+4*a*d+47*a*b*c+18*c^2*d+28*c*d^2+28*a*b^2*d+20*a^2*d^2+53*a*d*e^2+7*a^2*b^2*e+21*a*c^3*e+58*b*c*d^2*e+58*c^2*d^2*e+13*a^2*c*e^2+47*b*c*d*e^2+61*b*d*e^3,
46+70*b+76*c+3*a*c*d+50*c*d^2+2*a^2*b^2+34*b^3*d+34*b*c*d*e+2*d^2*e^2+21*b^4*c+44*a^3*c*d+54*a*b^2*c*d+39*b^3*c*d+27*a*b*c^2*d+6*a*c^3*d+83*c*d^4+66*d^5+55*b^4*e+66*b^2*c^2*e+27*a*b*c*d*e+43*c^2*d^2*e+45*a*b^2*e^2+52*b^2*d*e^2]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(7),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[22+50*b*c+74*a*d+7*d^2+33*b*e^2+37*d*e^2+73*b^2*c^2+33*b^2*d^2+47*c*e^3+59*a^3*b^2+30*a*c*d^3+62*c^2*d^3+68*a*d^4+35*a*b^2*e^2+7*b*c*e^3+89*a*e^4,
31*b*d+42*a*b*c+50*a*c^2+89*b*c*d+68*b^2*c^2+27*a*b*c*d+79*b*c^2*d+44*b*c*e^2+80*b*d*e^2+71*a*e^3+14*b^5+9*a^2*c^3+5*a*c^4+2*a^3*b*d+3*a^3*d^2+26*c^3*d^2+40*b*d^4+53*a^2*c^2*e+24*a^2*d^2*e+22*d^4*e+79*a*b*c*e^2+33*b*c^2*e^2+76*c^2*d*e^2,
33*b+50*e+61*a*b*c+74*c^3+52*a^2*e+8*a*b*e+75*a*d^3+19*a*b*d*e+7*a*d^2*e+36*d^3*e+55*a^2*e^2+57*a*c*e^2+83*c*d*e^2+16*b*e^3+13*a^3*b^2+33*b*c^4+4*c^5+51*a^2*b*d^2+82*b^2*c*d^2+62*a*b^2*c*e+47*c^3*d*e+90*a*b*d^2*e+55*a^2*d*e^2+14*a*b*e^3+71*d*e^4]))
    push!(
        ideals,
        Singular.Ideal(
            R,[85*b+55*c+18*b*c+66*b*d+70*a^3+7*a*b*c+79*d^2*e+33*a^2*b^2+31*a*b^3+55*a*b^2*c+36*d^4+80*a^3*e+29*b^2*c*e+74*d^3*e+12*c*e^3+78*a^2*b^2*c+25*a*b^2*c^2+53*c^3*d*e+15*a*c*d^2*e+3*b^3*e^2+69*a*c^2*e^2+67*c*e^4,
73+64*a*b^2+75*c^2*d+81*b^4+70*a^2*c^2+58*a^3*c^2+71*a*b^3*d+77*a^4*e+83*a^3*c*e+83*a*b^2*c*e+50*b*c*d^2*e+23*a^2*c*e^2+45*a*b*d*e^2+33*a*d*e^3,
11*a*b+24*a*e+a*c*d+29*b*d^2+70*b^2*e+33*c^2*e+26*b^3*c+54*b^2*c^2+73*a*b*c*d+59*b*c^2*d+11*d^4+39*d^2*e^2+42*c*e^3+15*a^2*b^3+46*a^4*c+7*a^3*b*c+51*b^3*c^2+60*a^4*d+64*a^2*c^2*d+64*a*b*c*d^2+85*b^2*d^3+30*d^5+23*a^4*e+83*b^3*c*e+41*a*b^2*d*e+86*b*c*d^2*e+2*a*b^2*e^2+23*b*c*e^3]))
    push!(
        ideals,
        Singular.Ideal(
            R,[3*a+65*a*d+7*a*c^2+44*a^2*e+74*c*d*e+32*e^3+33*c^4+89*b*c^2*d+49*b^2*d^2+37*c^2*d^2+46*a^2*e^2+59*c^2*e^2+22*b^3*c^2+72*a*b*c^3+5*a*b^3*d+61*b^4*d+82*a*b^2*c*d+33*d^5+49*b^2*c^2*e+71*c^4*e+11*a^3*d*e+60*b*c^2*d*e+78*a^2*d^2*e+42*b*c^2*e^2+24*b*d*e^3,
19*a*e+60*c^3+10*c^2*d+30*b*c*e+11*d^2*e+3*b*e^2+57*a^3*c+64*a*b*d^2+88*c*d^3+80*b*c^2*e+86*b*c*d*e+53*c*d^2*e+80*a*b*e^2+30*b*e^3+74*a^3*b^2+86*a^2*b^3+71*a^2*c^2*d+8*a*b*c^2*d+46*a*b*c*d^2+57*a^2*d^3+37*b*d^4+9*b*d^2*e^2,
44+3*a+53*b*d+5*b*e+57*a^3*d+a*b^2*d+46*a*c*d^2+65*a^3*b*c+12*b^2*c^3+64*b^2*c^2*d+82*c^4*d+53*a^2*c*d*e+30*b^2*d^2*e+65*a*d^3*e+40*b*c^2*e^2+6*a*c*d*e^2+53*e^5]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(8),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[80*b*d+35*c^2*d+39*b*c*e+11*c^2*e+49*a*d*e+60*a^2*b^2+67*a*b^2*c+66*a*c^2*d+11*c^3*d+4*c^2*d^2+43*a^3*e+52*b^2*c*e+55*c^2*d*e+48*b*d^2*e+43*c*e^3+61*b^3*c^2+64*a*c^3*d+24*a*b*c*d^2+10*a*c^2*d^2+61*b*d^4+5*c*d^4+7*a*b*d^2*e+60*a^2*b*e^2+90*b^3*e^2+18*c^2*d*e^2+23*a*d^2*e^2+38*b*d*e^3,
a+27*a*d+75*d^2+88*b*c^2+16*a^2*d+34*b*d^2+53*a*e^2+37*a^2*b^2+39*a*b*c^2+89*a^2*b*d+14*a^2*d^2+27*b^3*e+40*a*b^3*c+33*a*b^3*d+81*b^4*d+89*a*c^3*d+27*a^4*e+72*b^2*c^2*e+65*b^2*c*d*e+38*b^2*d*e^2+46*c*d^2*e^2+83*a*c*e^3+36*d*e^4,
45+72*d+50*a^2+3*c*e+15*a^3+10*b*e^2+35*a^3*b+5*a*b^2*d+59*a*b*c*e+6*a^2*e^2+83*a*b^3*c+18*b^4*d+9*a^3*c*d+63*b^3*d^2+51*a^2*c*d^2+4*b^2*c*d^2+11*c^3*d^2+18*a*b^3*e+89*a^2*c^2*e+13*a^2*c*d*e+40*b*d^3*e+87*c^2*e^3]))
    push!(
        ideals,
        Singular.Ideal(
            R,[16*a*c+31*b*d^2+24*a^2*e+70*a*d*e+15*d^2*e+71*d^4+b^2*c*e+58*a*c*e^2+33*b^4*c+74*a^3*c^2+86*a^2*b^2*d+7*c^3*d^2+53*b^2*d^3+86*a^4*e+48*a^2*b*c*e+64*a*b*c^2*e+88*b^2*c*d*e+42*a*c*d^2*e+a^2*b*e^2+4*a^2*d*e^2+6*c*d^2*e^2+60*b*d*e^3+7*e^5,
28*a+6*c*d+54*a^2*c+14*b^2*c+10*a^2*b^2+23*b^2*c^2+48*b^3*d+13*b*c^2*d+43*a*b*d^2+41*b^2*c*e+41*b*c^2*e+47*b^2*d*e+24*b*d^2*e+87*a^2*e^2+65*b^5+85*b^4*c+87*b^3*c*d+44*b^2*c^2*d+61*a^2*b*d^2+85*d^5+90*b*c^3*e+31*a^3*d*e+67*b*c^2*e^2+34*b*d^2*e^2+14*d^3*e^2,
46+52*d+20*a^2+56*b*e+9*d*e+88*a^3+24*a^2*b+25*a*d*e+19*b*d*e+20*a^3*b+58*b^4+40*b*c^2*d+66*a^2*d^2+14*b^2*d*e+67*a*e^3+26*a^3*b^2+55*a^3*c^2+68*a^2*b*d^2+16*a*b^3*e+4*c^2*d^2*e+5*b^2*c*e^2+89*c^3*e^2+30*a*b*e^3+31*c*e^4]))
    push!(
        ideals,
        Singular.Ideal(
            R,[31+89*b+9*a*c^2+76*b*e^2+51*a*b*d^2+52*c*d^3+64*a*b*c*e+77*a*b^4+32*a*b^3*d+43*b^3*c*d+62*d^5+4*a^2*b^2*e+62*a*b^2*e^2+86*a*c^2*e^2+36*a*b*d*e^2+65*a*c*e^3+72*c^2*e^3+90*a*d*e^3,
6*a^2+48*a*c+83*a*e+21*b*e+84*a*d^2+82*d^3+43*b^2*e+5*c*d*e+42*a^4+40*a^3*c+31*a^2*b*d+59*b^3*d+57*a^2*d^2+9*a*b*d*e+67*c*d^2*e+3*d^3*e+27*a*b*e^2+60*a^5+36*b^5+18*b^4*c+24*a^3*c^2+89*a^4*d+68*c^4*d+18*a*c*d^3+79*d^5+44*b^3*c*e+49*a*b*c^2*e+53*a^2*d*e^2+23*b*d^2*e^2+85*b*c*e^3+48*e^5,
87*a+39*b*e+3*a*b^2+3*b*d^2+44*c*d^2+84*b*e^2+57*b^3*c+72*a*b^2*d+5*a*d^3+17*a^2*b*e+25*a*b*c*e+11*c*d^2*e+70*d*e^3+71*a^5+58*a^4*c+35*a^3*b*d+44*a*b^2*c*d+a^3*b*e+20*a*b^3*e+a^2*c*d*e+64*a*c^2*d*e+41*a^2*b*e^2+2*c*d*e^3]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(9),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
    push!(
        ideals,
        Singular.Ideal(
            R,[53+6*a+4*a*d+15*c*d+5*d^2+20*a*e+37*c^2*d+4*d^3+43*b^2*e+55*d*e^2+55*c^4+21*b^2*d^2+11*a*b*c*e+67*a*c*e^2+12*a*e^3+47*b^2*c^2*d+35*b*c^3*d+25*b^4*e+80*a^3*c*e+8*b^2*c*d*e+45*a^2*d^2*e+11*b*c*d^2*e+41*a*b*c*e^2+73*c^2*e^3+55*b*d*e^3,
45*a*c+44*b*c*e+77*a*d*e+21*c*d*e+16*a*e^2+42*d*e^2+8*a^3*c+5*a*d^3+10*b*c*d*e+2*c*d^2*e+66*d^3*e+64*c^2*e^2+19*b^5+8*a^2*b^2*c+71*a^4*d+58*b^4*d+39*a*b^2*c*d+25*b^2*c^2*d+41*a*b*c*d^2+40*b^2*d^3+79*c^2*d^3+49*a^2*b^2*e+77*a*b^2*c*e+26*a*b*c^2*e+64*b^2*c*d*e+30*a*c*d*e^2+87*e^5,
42*b+48*a^2*e+64*b^2*e+28*a^3*b+72*b^3*c+64*b*c^3+24*c^4+77*a^3*d+17*b^3*d+18*a*c^2*d+69*a*c*d^2+41*b^2*c*e+8*c*d^2*e+71*a^5+51*a^2*b^3+8*b^5+54*a^2*c^3+24*b^2*c^2*d+9*c^4*d+25*a^3*d^2+50*b^3*d^2+24*a*b^3*e+88*b^2*c^2*e+41*b*c^3*e+48*a*b*d*e^2+40*c*d^2*e^2]))
    push!(
        ideals,
        Singular.Ideal(
            R,[14*a^3+9*b^2*c+84*a*b*d+17*a*d^2+67*c^2*e+19*b*d*e+60*a^3*c+17*a*c^2*d+10*b*c^2*d+49*a*d^3+69*a^3*e+6*b^2*c*e+67*a*c*d*e+40*b*c*d*e+74*a^3*b*c+37*a*b^3*c+55*b^3*c^2+30*a^2*b^2*d+16*a^2*b*c*d+18*b^3*c*d+85*a^2*c^2*d+90*c^4*d+21*a^2*b*d^2+11*a^2*c*d^2+59*b*c^3*e+49*a*b^2*d*e+42*a^2*c*d*e+13*a*c^2*d*e+41*b^2*d^2*e+60*b^2*d*e^2+49*b*d*e^3,
34+12*b+42*a^2+3*b^2+15*a*d+10*a*b*c+38*a*b*d+68*a*c*d+82*a^2*e+39*b*e^2+16*c^4+90*a*c*d^2+26*d^4+41*a^2*c*e+11*b^3*c^2+69*a^2*b*c*d+85*c^4*d+19*b*c^2*d^2+73*a^3*b*e+19*a^2*b*c*e+87*a*b*d^2*e+4*a*b*d*e^2+25*b^2*d*e^2+35*c^2*e^3+15*b*e^4,
79*b+87*b*c+27*c*e+78*a^2*c^2+32*c^3*d+7*b*c*d^2+43*c*d^3+30*a*b*c*e+31*b*c^2*e+29*c^3*e+57*a^2*d*e+78*d^3*e+12*a*b^4+49*a*b^3*c+28*a*c^4+80*b*d^4+42*b^4*e+55*a^2*b*c*e+32*b^3*d*e+87*b*c^2*d*e+42*a^2*d^2*e+80*b*c*d^2*e]))
    push!(
        ideals,
        Singular.Ideal(
            R,[6*c+22*b*c+11*a^3+13*b^2*c+84*a^2*d+32*a*b*d+2*a*d*e+15*a^2*b^2+49*b*c^2*d+61*b^3*e+40*a*c^2*e+12*b^2*d*e+22*c*d^2*e+72*b^2*e^2+26*b*c*e^2+54*a^3*c^2+87*a^2*b*c*d+46*c^4*d+17*a^2*d^3+17*b^3*c*e+66*b^2*c^2*e+60*b*d^3*e+12*a^2*b*e^2+32*c^2*e^3+12*b*d*e^3+57*c*d*e^3,
72+81*b+45*b*c+67*e^2+38*b*d*e+15*a*e^2+36*c*e^2+21*b^3*e+75*a*c*d*e+8*a*d^2*e+57*b*d^2*e+77*d^3*e+4*a^4*b+10*a^3*b^2+56*b^5+76*a^2*b*c^2+28*a*c^4+53*a^2*b^2*d+37*b^3*d^2+71*a*b*c*d^2+17*a^2*b^2*e+84*b^4*e+63*a*b*c^2*e+42*b^2*c*d*e+49*a*d^3*e+77*a*b^2*e^2+66*a^2*c*e^2+72*b^2*d*e^2+64*c*d^2*e^2+12*c*d*e^3,
24*d^2+85*e^2+80*b^2*d+9*c^2*d+42*a*d*e+48*a^3*c+75*a^2*b*d+81*c^3*d+6*a*b*d^2+23*c*d^2*e+47*d^3*e+53*b*c*e^2+36*e^4+10*a*b^2*c^2+61*a*b^3*d+17*a^3*c*d+46*a*c^2*d^2+5*a*b^2*d*e+65*c*d^2*e^2+89*b^2*e^3+48*d^2*e^3+42*c*e^4]))
    for i in ideals
        runAll("sparseid(3,0,5,100-(1),90)", i, S, StartOrd, TarOrd)
    end
    ideals =[]
end
