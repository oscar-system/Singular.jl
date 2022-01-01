include("GroebnerWalkFinalBenchmarkProcedures.jl")
include("runbenchmark.jl")
include("GroebnerWalkFinalBenchmark.jl")
include("runbenchmark2.jl")
include("BenchmarkHelper")
include("readWriteHelper.jl")
include("Examples")
using DataFrames
using CSV

function benchmarkVar4()
    cd("/Users/JordiWelp/Results")
    prepare()
    prepare2()
    prepareAlloc()
    dim = 4
    ve = [1, 1, 1, 1]
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (a, b, c, d) = Singular.PolynomialRing(
        Singular.QQ,
        ["a", "b", "c", "d"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)
    ideals = []
            push!(
                ideals,
                Singular.Ideal(
                    R,87+55*b+17*a^3+82*a^4+17*b^2*c^2*d+71*a*b^2*d^2,
        26*c*d+14*a^3+89*b^4*c,
        88*a^2*b*c+61*a^2*d^2+48*c^5))
            push!(
                ideals,
                Singular.Ideal(
                    R,27*a^2*b^2*d+60*a*b*c*d^2,
        23+31*b*d^2+42*a^2*c*d+32*b*c^2*d+29*a*b^3*c,
        d+6*a*d+67*a^3+48*a*b^3+86*a^2*b*d^2))
            push!(
                ideals,
                Singular.Ideal(
                    R,85*c*d^2+10*d^3+11*a^3*c+31*a^2*c*d+61*a^2*c*d^2,
        57+86*a^2*d^2+86*b^3*c^2+58*b^2*c^2*d,
        76*d+90*c*d+40*a*c^3*d))
            for i in ideals
                runAll("sparseid(3,0,5,100-(1*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,3+35*a*d+78*c*d+47*a*c^2+30*c*d^2+65*a*c^3+55*a*b*c*d+25*a*b^3*c+52*b^2*c^3+52*a*c*d^3,
        71*a^3+65*a^4+9*a^2*b*c+65*b*c^2*d+25*a^4*c+32*a^3*c*d+45*a^2*c*d^2,
        49*d+16*a^2*b^3))
            push!(
                ideals,
                Singular.Ideal(
                    R,38*b^4+80*a^2*b^2*c+41*b*c^3*d+7*b^2*d^3,
        4+57*d+44*b*c+50*a^4+24*a*b^3+29*a^3*c+55*b*c^3*d+82*a*c^2*d^2+77*c^3*d^2,
        79*d^2+13*a*c^2+82*c^3+14*a*c*d+37*b^2*d^2+68*b^2*d^3))
            push!(
                ideals,
                Singular.Ideal(
                    R,70*b+4*a^2*c+85*c^4+21*a^2*b^2*c+17*a*b^3*c+31*a^2*c*d^2+14*c^3*d^2,
        2*c^2+35*d^2+44*b^3+34*a*b^2*c+87*b^2*c*d+31*a*c^2*d+6*a*c^4,
        16+88*a*c^2+50*c^2*d^2+70*a^3*b*c+62*a^3*c^2))
            for i in ideals
                runAll("sparseid(3,0,5,100-(2*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,67*a*c*d+3*b*d^2+85*c*d^2+71*a^2*c^2+45*a*d^3+39*b*d^3+77*b^5+42*a^2*b*c*d+54*b^3*c*d+3*b^2*c*d^2+35*b*c^2*d^2+65*a*c*d^3,
        40+47*b+66*a*b+47*a*b^2*d+78*a*b*c*d+41*c^2*d^2+80*a^4*d+7*c*d^4,
        49*d^2+89*a^3+36*a*b*c^2+24*a^2*c^3+82*a*b^3*d+33*b^3*d^2))
            push!(
                ideals,
                Singular.Ideal(
                    R,62*d^2+44*a^3+43*a*b*d^2+13*a*c*d^2+14*b*c*d^2+65*a^3*b*d+74*a^2*c^2*d,
        5+57*b+14*a^3+65*a^2*b+79*a*c^2+65*a^3*d+82*b^2*c*d+44*a^2*b^2*c+55*b^4*c+34*a*b^2*d^2,
        30*b^2+29*b*c^2*d+49*b*c*d^2+22*a^3*b*c+21*a^2*c^3+19*a^2*b^2*d+61*b^2*c*d^2+50*a*b*d^3+9*b*d^4))
            push!(
                ideals,
                Singular.Ideal(
                    R,46+79*b^2*c^2+54*a^5+5*a*b^4+34*a*b^2*c^2+19*a^2*b*d^2+11*a*b^2*d^2,
        69*b*c+86*a^3+52*c*d^2+87*b^2*c^2+82*a*c^2*d+51*a^2*b^2*d+72*b^2*c^2*d+42*d^5,
        77*b+77*d^2+50*a^2*c+55*a*c^2+67*b^4+b^2*c^2+17*a*c^3+39*b*c^2*d+70*a^3*b*c+50*b^3*c^2+24*a*b*c^3))
            for i in ideals
                runAll("sparseid(3,0,5,100-(3*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,6*c^2+49*b*d+54*a^3+35*a*b^2+6*c*d^3+28*c^5+7*a*c^3*d+3*b*c*d^3,
        80*a^2*b+82*a*b^3+8*b^3*d+6*b^2*d^2+42*b*c*d^2+80*a^5+42*a*b^4+50*a^4*d+69*a*b*c^2*d+87*c^3*d^2+13*a*c*d^3+72*b*d^4,
        27+16*b+81*c*d+53*a^3+41*b^2*c+63*a*b^2*c+21*a^3*d+5*c^3*d+12*a^2*d^2+36*a^3*b^2+85*b^2*c^3+28*a^2*b*c*d+77*b^2*c^2*d))
            push!(
                ideals,
                Singular.Ideal(
                    R,33*a+6*c*d+38*a*b^3+82*a^3*c+61*a^2*b*c+48*b*c^3+32*a^4*c+25*a^2*b^2*c+34*a^3*b*d+54*a*b^3*d+90*a^2*b*d^2+74*a*c^2*d^2+3*c*d^4,
        55+80*b^2+74*d^2+77*a*b^2+71*b^2*c+63*c^3+b*c*d+75*b^4+32*a^3*d+46*c^2*d^2+73*a*b^2*c^2+40*c^5+43*a^2*b^2*d+43*b^4*d+75*a^2*d^3,
        51*a^2*d+10*b^3*d+52*c*d^3+61*b^4*c+77*a^2*b*c^2))
            push!(
                ideals,
                Singular.Ideal(
                    R,80*c+68*b^2+52*b*c^2+53*b*d^2+14*d^3+63*b*c^3+34*a^2*c*d+42*b^2*c^2*d+14*a*c*d^3,
        11*b^2+75*a*b^2+77*a^2*b^2+88*b^3*c+71*b*c^2*d+74*a*b^2*c*d+51*a*c^3*d+67*b^2*d^3+68*a*d^4+83*b*d^4,
        38+28*b*d+22*b*d^2+27*a*b^3+39*b^2*c^2+26*a*b^2*d+81*c^3*d+16*a^5+81*b^3*c^2+45*a*c^4+27*a^4*d+80*a*b^3*d+11*a*c^2*d^2+21*b*c*d^3))
            for i in ideals
                runAll("sparseid(3,0,5,100-(4*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,25*a^2+7*a*b^2+5*a^2*c+25*c^3+83*a*b^2*d+49*a^2*c*d+72*a*b*c*d+32*c*d^3+30*a*c^4+53*b^3*d^2+37*a^2*c*d^2+73*b^2*c*d^2,
        67*b^2+48*a*c*d+55*b*c*d+61*a^2*c*d+24*a^3*b*c+32*a^2*c^3+78*a^4*d+33*a*b*c^2*d+29*b^3*d^2+40*a*b*c*d^2+36*c^3*d^2,
        36+72*a+11*d+77*a*d+85*b^2*d+4*a^3*b+65*a^2*c*d+15*b^2*c*d+70*b*c^2*d+52*c^3*d+57*c*d^3+38*a*b^4+6*a^2*b^2*c+15*a*b*c^2*d+48*c^4*d+83*a*b^2*d^2+78*b^2*c*d^2))
            push!(
                ideals,
                Singular.Ideal(
                    R,5*c*d+48*a^3*b+29*a*b^3+30*a^2*c^2+7*a*b*c^2+29*b*c*d^2+80*a*b^4+75*b*c^4+35*b^2*c*d^2+14*a^2*d^3+51*b*d^4+82*c*d^4,
        23+18*b+64*d+70*c*d+67*d^2+12*a^2*b+36*a*c*d+19*b^3*c+50*b*c^2*d+73*a*c*d^2+42*a^4*c+7*a^2*c^3+9*a^2*c^2*d+23*a*b*c^2*d+88*b^2*c^2*d+68*c^4*d,
        22*a^3+40*a^2*b+59*a*b^2+10*a*b*d+90*a*b^2*c+36*a*c^3+33*a*c^2*d+75*a^2*b^3+66*a^3*b*c+60*c^3*d^2+29*a*b*d^3+2*b*d^4))
            push!(
                ideals,
                Singular.Ideal(
                    R,44*d+42*b*d+64*d^2+40*a*b*c+27*a^3*b+27*a*c^3+73*a*b^2*d+69*a*c^2*d+55*b*c*d^2+58*a*d^3+74*a^4*b+35*a*b*c^3+3*b^3*c*d+68*a*c^3*d+76*c^4*d+66*c^2*d^3,
        56*a^2*c+83*a*b*c+17*a*c^2+61*c^2*d+40*b^3*c+a*b*d^2+2*b^2*c^3+43*a^4*d+65*b^3*c*d+77*a*b*c^2*d+45*a^3*d^2+79*b^2*c*d^2+12*c*d^4,
        26+24*c+62*a^2+70*a^2*b+4*b^4+32*b*c^3+75*c^2*d^2+29*b^4*d+25*b^3*c*d+79*b*c^3*d+75*b^2*d^3))
            for i in ideals
                runAll("sparseid(3,0,5,100-(5*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,17*c^2+24*b^3+71*c*d^2+38*a^3*b+90*a*b^3+90*a^2*b*c+40*a*b^2*c+27*b^2*c*d+70*a*c^2*d+52*a^3*b*d+32*a^3*c*d+12*a*b^2*c*d+78*b*c^3*d+86*c^4*d+7*a*b^2*d^2+50*a*d^4,
        44*a*c^2+2*a*b*d+80*c*d^2+60*d^3+87*a^4+85*a^2*b*c+88*a*c^3+67*b^2*d^2+17*b^5+68*a^2*b*c^2+64*b^2*c^3+86*b^2*d^3+49*c^2*d^3+46*b*d^4,
        25+77*c+77*d+75*a*d+52*c*d+89*d^2+11*b^3+65*a*c*d+64*b^3*c+10*a*b*c^2+22*b^2*d^2+14*a^5+60*a^3*b^2+22*a^2*b*c*d+34*a*b^2*c*d+81*a*b^2*d^2+20*b^3*d^2+7*b*c^2*d^2+84*a*c*d^3))
            push!(
                ideals,
                Singular.Ideal(
                    R,41*a^2+12*a*b^2+79*b^3+25*a^2*d+44*c^2*d+39*a^3*b+77*a*c^3+62*b*c^3+72*b*c^2*d+66*b^2*c^3+45*b^3*c*d+3*b*c^3*d+24*c^4*d+64*b^3*d^2+37*b*c^2*d^2+76*c*d^4,
        17*a^2+12*b^2+30*a^2*b+72*b*c*d+81*a*d^2+31*a^4+31*b^4+82*a^3*c+26*a^2*c^2+84*a^3*d+62*a^2*b*d+50*a*b^4+21*b^5+15*a^2*b^2*c+6*a*b^2*c^2+84*a*b*c^3+86*a*b^2*c*d+56*b*c^2*d^2+64*a*b*d^3))
        88+90*a+88*b+30*c*d+19*b*d^2+86*a^2*b*d+88*b^3*d+53*b^2*c*d+62*a^3*b*c+41*b^3*c^2+28*a*b^2*d^2+48*b^3*d^2+83*a*b*c*d^2+30*c^3*d^2,
            push!(
                ideals,
                Singular.Ideal(
                    R,23*a*d+45*a*b^2+9*a*d^2+53*d^3+85*b^2*c*d+42*a*c^2*d+25*a^2*d^2+12*a*b*d^2+29*b^2*d^2+45*b*c^4+48*a*b^2*d^2+52*b^3*d^2+36*a^2*d^3+14*b*d^4,
        52*c+6*c^2+24*c*d+42*a^2*b+24*a^4+41*a*c^2*d+16*a^3*b*c+55*a^2*c^3+44*b^2*c^3+66*a^2*b^2*d+11*a^3*c*d+10*a*c^2*d^2+4*a*c*d^3+87*b*c*d^3,
        59+2*d+74*c*d+b^3+83*a*b*c+73*a*b*d+a*c*d+69*a^2*b*c+40*c^4+30*a^3*d+58*a^2*c*d+14*c^3*d+22*a*d^3+35*a^3*b^2+79*b^2*c^3+15*a^3*c*d+68*a*c^3*d+10*c^4*d+31*a^2*c*d^2+44*b^2*c*d^2+17*c^2*d^3))
            for i in ideals
                runAll("sparseid(3,0,5,100-(6*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,2*c^2+25*b^2*c+79*b*c^2+60*a^3*b+47*b^2*c*d+28*b*c*d^2+28*c^2*d^2+38*a*b^4+82*a^2*b*c^2+54*a^4*d+73*a*b^2*d^2+77*a*c^2*d^2+51*d^5,
        1+33*d+41*a*b+79*c*d+9*a^3+41*a*b*d+43*b*c*d+51*a^4+59*b^2*c^2+72*a*b*c*d+35*a^2*d^2+70*a*b*d^2+63*b^2*d^2+11*a^2*b^2*d+54*b^4*d+50*b^3*c*d+45*a^2*c^2*d+28*b*c^3*d+14*a^2*b*d^2+15*c^3*d^2+84*a*b*d^3+11*c^2*d^3,
        78*a+65*a*b+78*a*d+29*a*b^2+30*a^2*c+43*b*c^2+71*a*b*d+36*a^2*b*c+85*a^3*d+51*b^3*d+49*a^2*c*d+75*d^4+17*a^4*b+70*a*b^4+4*a^3*b*c+59*b*c^4+32*a^3*c*d+2*c^4*d+61*b^3*d^2+45*a*b*d^3+60*b*c*d^3))
            push!(
                ideals,
                Singular.Ideal(
                    R,40*b^2+53*c^2+13*b*c*d+54*a^4+84*a^2*c^2+83*b*c^3+38*a*b*d^2+18*a^2*b^3+76*b^3*d^2+23*b*c^2*d^2,
        10*d+63*b*c+67*b^3+69*a^2*c+84*a*c^2+31*c^3+88*b^2*d+66*b*d^2+31*b^4+42*a^3*d+46*b*c^2*d+14*c^3*d+64*a^5+30*a*b*c^3+56*b^2*c^3+37*a^2*b^2*d+85*a^2*c^2*d+10*a*c^3*d+82*c^4*d+90*a^2*b*d^2+2*a*b*d^3+87*a*c*d^3+76*b*d^4,
        50+64*d+37*a*b+9*b*c+8*a*b^2+19*c^3+68*a^4+58*a^2*c^2+68*a^3*d+34*a*b*c*d+87*b^2*c*d+64*a*c^2*d+22*a*d^3+76*a^3*b*c+82*b^4*c+53*a^2*b*c^2+10*a*c^4+81*b*c^4+57*c^5+49*a^4*d+25*a^2*b*d^2+79*b^2*d^3+47*c^2*d^3))
            push!(
                ideals,
                Singular.Ideal(
                    R,59*a+21*a*c^2+65*d^3+42*a*b^3+24*a*b*c^2+71*a*c^3+13*a*b*c*d+86*b^2*d^2+86*b^3*c^2+70*b^2*c^3+12*a*b^3*d+68*a^3*c*d+21*b^3*c*d+27*a^2*b*d^2+77*a*b*d^3+49*c^2*d^3+25*b*d^4+43*c*d^4,
        68+20*a*b+11*a*d+71*b^3+36*a*d^2+49*a^3*b+78*b^4+88*a^3*c+49*a*b^2*c+16*c^2*d^2+8*a^3*c^2+61*a^3*b*d+25*b^2*c^2*d+25*b^3*d^2+21*a^2*c*d^2+36*b^2*c*d^2+78*a*c^2*d^2,
        12*b+75*c^2+68*c*d+46*d^2+68*a^2*b+70*a*b^2+59*b^2*c+33*a*b*d+83*c^2*d+7*a^2*b*c+58*a*b*c*d+70*c^3*d+46*b*c*d^2+12*c*d^3+68*a^3*c^2+49*a^4*d+37*b^2*c^2*d+a^2*c*d^2+10*a^2*d^3+7*a*b*d^3+16*b*c*d^3))
            for i in ideals
                runAll("sparseid(3,0,5,100-(7*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,89*c+37*c^3+29*b*c^3+70*a^2*c*d+38*a*c^2*d+63*a^2*d^2+73*b*d^3+40*c*d^3+59*d^4+69*a^2*b^3+86*b^4*d+26*a*c*d^3+14*c*d^4+49*d^5,
        31*d+50*a*c+56*c^2+28*c*d+46*a*b*d+74*a*c*d+12*a*d^2+57*d^3+86*a^2*b*c+6*b^2*d^2+72*a*c*d^2+26*b*d^3+20*a^4*b+83*a^3*b*c+27*b^3*c^2+a*b*c^3+67*b^2*c^3+64*a^4*d+19*a^2*c^2*d+71*a^2*c*d^2+28*a*b*c*d^2+89*b*c^2*d^2+56*c^2*d^3+74*c*d^4,
        33+6*b*c+58*d^2+87*b^2*c+3*a*b*d+11*a*c*d+64*a*d^2+83*d^3+28*a^4+11*a*b^3+21*a*b^2*c+69*b^2*c*d+44*c^2*d^2+33*a*d^3+4*b^4*c+50*a^3*c^2+29*b^3*c^2+60*a^3*b*d+63*a^2*b^2*d+11*a*b^2*c*d+82*a^2*b*d^2+43*b^3*d^2+89*c^3*d^2+16*d^5))
            push!(
                ideals,
                Singular.Ideal(
                    R,77+57*d+17*d^2+69*a^2*b+67*a*c^2+43*a^2*d+23*a^2*b^2+45*a^2*c^2+32*c^4+50*a*b^2*d+75*b^2*d^2+88*a^5+22*b^5+66*a^2*b*c^2+68*b*c^4+7*a^3*b*d+64*a^3*d^2+9*b*c^2*d^2,
        74*a^2+16*b*c+45*b^2*d+84*c^2*d+61*d^3+69*a^2*c^2+14*a*c^3+89*b*c^3+74*a*b*c*d+83*a*c^2*d+80*a*b*d^2+38*a^4*b+51*b^3*c^2+33*b^2*c^3+54*c^5+86*a^3*b*d+69*b^2*c^2*d+5*c^4*d+89*a^2*b*d^2+50*b^2*c*d^2+42*b*c^2*d^2+35*c^3*d^2+16*b^2*d^3,
        20*b+46*a*b+7*c^2+33*a^2*b+56*a^2*c+34*b^2*c+12*b*d^2+16*a*b^2*c+35*a^2*c^2+85*b*c^3+16*a^2*b*d+17*b*c*d^2+74*b*d^3+4*a*b^3*c+47*a*b*c^3+34*c^5+84*a^3*b*d+64*a*b^2*c*d+a^2*b*d^2+62*c^3*d^2+18*b*d^4))
            push!(
                ideals,
                Singular.Ideal(
                    R,54+47*b+8*a*c+30*c^2+6*a*b^2+16*b*d^2+56*c*d^2+77*a^3*c+2*a*b^2*c+71*a^2*c^2+83*a*b*c^2+13*b*c^3+63*b*c^2*d+23*b*c*d^2+68*a*d^3+65*b*d^3+82*a^5+7*a^3*b*c+58*b^4*c+25*a^2*b*c^2+81*b^2*c*d^2+65*b*c^2*d^2+12*c^3*d^2+76*b*c*d^3,
        80*a*c^2+87*a*b*d+72*c*d^2+48*b^4+72*b^3*c+75*a*b*c^2+89*a^5+25*b^5+31*a^4*c+48*a^3*b*c+55*a^2*b^2*c+56*a*b^3*c+84*a^2*b*c^2+54*a^3*d^2+46*b^2*c*d^2+19*b*c*d^3,
        62*a+66*a*b+53*a*c+2*c^2+64*a*b^2+19*b^2*c+48*c^3+42*c^2*d+77*b^2*c^2+39*a^2*c*d+50*a^2*d^2+83*b*c*d^2+14*a*d^3+26*a^5+90*a^4*b+28*a^4*c+42*b^2*c^3+28*b*c^4+25*c^5+64*a*b*c^2*d+5*a^2*c*d^2+62*b*c^2*d^2))
            for i in ideals
                runAll("sparseid(3,0,5,100-(8*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,49*d+54*a^2+68*c^2+19*a^2*c+23*a*b*c+31*a*c^2+81*a*c*d+70*b*d^2+28*c*d^2+44*a^3*c+78*b^2*c^2+87*a*c^3+46*c^4+20*b*c^2*d+34*a*b*d^2+32*c^2*d^2+32*a^5+52*c^5+27*a^4*d+48*b^4*d+36*b^3*c*d+16*c^4*d+47*c*d^4,
        11*a+79*c*d+14*a^2*b+27*a*b*c+57*a*d^2+90*a^3*b+43*b^3*c+71*a^3*d+70*a^2*c*d+88*b^2*c*d+2*a*b^4+75*a^3*c^2+14*a^2*b*c^2+57*b^3*c^2+79*a*c^4+28*b*c^4+79*a^3*b*d+83*a*b^3*d+44*a*b*c^2*d+50*b^2*c^2*d+10*c^2*d^3+24*b*d^4,
        29+52*b+5*a^2+86*a*b+83*b*d+83*a*b^2+75*a*b*d+80*a^3*b+6*a^2*b^2+70*a^3*c+72*b*c^3+60*a^2*c*d+16*b*d^3+63*d^4+56*a^4*b+74*a^4*c+84*a*b^3*c+80*b^4*c+81*a^3*c^2+87*a^2*b*c^2+87*a*b^3*d+49*b^3*c*d+81*a*c^3*d+59*b*c^3*d+71*a*d^4+78*c*d^4))
            push!(
                ideals,
                Singular.Ideal(
                    R,83*a+23*b*c+54*c*d+4*d^2+69*a^2*b+31*a*b*c+6*a^2*d+47*a*b*d+76*b^2*d+5*c^2*d+56*a*b^3+18*b^2*c^2+a^2*b*d+30*a*b*d^2+34*a*c*d^2+54*d^4+63*a^2*b^3+40*a^2*b^2*c+33*a^2*b*c^2+78*c^5+52*a^4*d+90*a^3*b*d+63*a^2*b*c*d+50*a*b*c^2*d+81*a*c^3*d+60*b^2*c*d^2+45*a*b*d^3+82*b^2*d^3,
        88*a+84*a^2+61*a^2*d+36*b*c*d+37*a^2*b^2+32*a*b^2*c+19*c^4+47*a^2*b*d+66*a^2*c*d+89*b^2*c*d+21*c^2*d^2+13*a*d^3+65*a^5+87*a^3*b^2+36*a^2*b^3+85*a*b^4+85*a*b^2*c^2+33*b^3*c^2+23*b*c^4+46*a^3*b*d+86*a^3*c*d+49*a*c^3*d+25*c^4*d+15*a^2*b*d^2+56*a*b*c*d^2+70*a*d^4,
        90+52*c+12*a*c+10*b*c+75*b^3+67*b*d^2+32*c*d^2+21*b^4+44*a*b^2*d+53*b^2*c*d+54*a*b*d^2+13*a*d^3+40*a^3*b^2+89*a^2*b*c^2+62*a^2*c^3+5*a^3*c*d+28*a*c^3*d))
            push!(
                ideals,
                Singular.Ideal(
                    R,63*a*c+b^3+27*a*b*c+6*b*c^2+76*c^3+89*b^4+36*a^2*b*d+63*b*d^3+70*a^5+80*b^5+45*b^4*c+a^2*b*c^2+14*b^3*c^2+5*a*b*c^3+3*c^5+54*a*b^2*c*d+17*c^4*d+60*c^3*d^2+69*a^2*d^3+87*c^2*d^3,
        43*a+23*d+23*c^2+43*d^2+71*a^2*b+35*b^3+75*a*c^2+12*a^2*b*c+12*a*b^2*c+2*b^3*c+63*b*c^3+11*a*b^2*d+90*b^3*d+28*a*b*c*d+20*b*c^2*d+34*c*d^3+90*a^5+26*a*b^4+66*b^5+2*a^4*c+33*a^2*b^2*c+27*a*c^4+42*a^2*b*c*d+b^3*c*d+65*a*b*c^2*d+86*a*c^3*d+87*b^2*c*d^2+60*c*d^4,
        72+55*c+76*b^2+30*a*c+85*b*d+14*a*b^2+34*a^2*c+12*b^2*c+57*a*c*d+25*a^3*b+65*a^2*b^2+29*a^2*b*c+50*a*b^2*c+32*c^4+46*a^2*c*d+53*a*b*d^2+72*a*b^3*c+43*a^3*d^2+78*a^2*c*d^2+78*a^2*d^3+81*b^2*d^3+52*b*d^4+73*c*d^4))
            for i in ideals
                runAll("sparseid(3,0,5,100-(9*2),90)", i, S, StartOrd, TarOrd)
            end
            ideals =[]
            push!(
                ideals,
                Singular.Ideal(
                    R,52*c+52*d+52*a^2+89*b*c+11*a^2*c+81*a^2*d+68*a*c*d+51*a^4+19*a^2*b^2+13*b^4+78*a^3*c+11*b^2*c^2+87*b*c^3+47*a*b^2*d+81*a*d^3+21*a^3*b^2+88*b^4*c+44*a^3*c^2+75*a^3*c*d+32*a^2*b*c*d+15*b*c^3*d+79*a^2*b*d^2+77*a*c^2*d^2+24*c^3*d^2+51*b*d^4,
        15*b+17*b*c+73*a^3+19*a*b*c+30*b^2*c+16*b*c^2+53*a*b*d+45*d^3+40*a^3*c+25*a*b^2*c+35*a*c*d^2+55*c^2*d^2+40*a^3*b*c+34*a*b^3*c+66*a^2*b*c^2+30*a*b*c^3+32*c^5+30*a*b^3*d+82*b^3*c*d+76*b^2*c^2*d+88*c^4*d+72*a^2*b*d^2+76*c^3*d^2+49*b*d^4,
        50+18*a*c+87*b*c+76*c*d+59*a^3+43*b^2*d+7*b*c*d+61*a^2*b^2+39*a*b^2*c+54*a*b*c^2+30*a*c^3+9*a^2*b*d+40*b^3*d+14*a^2*c*d+88*b*d^3+13*c*d^3+48*a^2*b^3+40*a^4*c+76*b^4*c+81*a*c^4+13*a^4*d+45*b^4*d+7*a^2*b*c*d+17*b*c^2*d^2+11*a*b*d^3+59*b*c*d^3+38*b*d^4+89*d^5))
            push!(
                ideals,
                Singular.Ideal(
                    R,29*b+68*a*d+6*b^3+27*a^2*c+28*b^2*c+a*d^2+27*d^3+23*a^3*b+37*b^3*c+28*a^3*d+43*b^3*d+7*a*b*c*d+59*c^3*d+44*a*b*d^2+67*a^2*c^3+62*a^2*b*c*d+75*a*b^2*c*d+52*a*c^3*d+60*a^2*b*d^2+55*b^3*d^2+72*a*b*c*d^2+69*a^2*d^3+4*c^2*d^3,
        11*c+12*c^2+38*a*b*c+68*a^2*d+3*b^2*d+54*b*d^2+80*a^2*b^2+69*a*b^3+a^2*b*d+35*b^3*d+87*a*c^2*d+80*a^2*d^2+6*b*d^3+46*a^2*b^3+23*a^4*c+84*b^3*c^2+50*b^2*c^3+33*a*c^4+5*a^2*b^2*d+69*a*b^2*c*d+89*a^3*d^2+28*a^2*b*d^2+34*b^3*d^2+22*a^2*c*d^2+52*a*b*c*d^2,
        7+35*b+69*b^2+4*a*d+74*b*d+21*d^2+39*a^2*d+9*a*b*d+53*b*d^2+66*a^3*b+86*a^2*b^2+49*a*b^3+58*a^2*b*d+51*b^3*d+87*a*c^2*d+22*a^2*d^2+27*a*b^3*c+18*a^2*b*c^2+45*a^2*c^3+12*a*b*c^3+4*a^3*b*d+78*a^2*b^2*d+19*b^3*c*d+23*a^2*c^2*d+56*c^4*d+26*b^2*c*d^2+13*a^2*d^3+87*a*b*d^3+65*c*d^4))
            push!(
                ideals,
                Singular.Ideal(
                    R,85*b^2+5*b*d+27*d^2+56*c*d^2+81*a^2*b^2+6*a*b^3+54*a^2*c^2+26*a*b*c^2+27*b^2*c^2+37*b*c^3+22*a*b^2*d+69*b*c^2*d+78*b^4*c+82*c^5+39*a^2*b^2*d+73*a*b^2*c*d+50*b^3*c*d+83*b^2*c^2*d+14*a*c^3*d+34*c^4*d+9*a^3*d^2+68*a*b^2*d^2+a*b*c*d^2+74*a^2*d^3+35*a*b*d^3+8*c*d^4,
        87*d+10*b^2+18*a*b*c+71*c^3+63*a^2*d+30*a*b*d+7*b^2*d+9*c^2*d+89*b*d^2+27*c*d^2+10*a*b^2*d+41*a^2*c*d+85*a*c^2*d+51*a^2*d^2+23*a*b*d^2+88*b^2*d^2+26*d^4+6*a^2*b^2*c+90*a*b^3*c+81*a*b^2*c^2+43*a^4*d+40*a^3*b*d+69*b^4*d+22*a*c^3*d+78*a^2*b*d^2+70*a*b^2*d^2+68*b^2*d^3+61*b*d^4+43*d^5,
        2+39*b+39*c+29*a*b+82*a*d+38*a^2*c+67*a^2*d+51*a*d^2+73*b^4+9*a^3*c+70*a*b*c*d+63*b^2*c*d+73*a*d^3+3*b*d^3+46*a^5+42*a*b^4+7*a^2*c^2*d+88*a^3*d^2+62*a^2*b*d^2+38*b^3*d^2+60*a^2*d^3+43*a*d^4))
        for i in ideals
            runAll("sparseid(3,0,5,100-(10*2),90)", i, S, StartOrd, TarOrd)
        end

end
