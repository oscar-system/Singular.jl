var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Getting-Started-1",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "section",
    "text": "Singular.jl is a Julia interface to the Singular computer algebra system. It was written by Oleksandr Motsak, William Hart and other contributors, and is maintained by William Hart, Hans Schoenemann and Andreas Steenpas. It is part of the Oscar project.https://www.singular.uni-kl.de/ (Singular website)\nhttps://github.com/wbhart/Singular.jl (Singular.jl source code)\nhttp://wbhart.github.io/Singular.jl/ (Singular.jl online documentation)The features of Singular so far include:Singular integers, rationals Z/nZ, Z/pZ, Galois fields\nMultivariate polynomials\nIdeals over polynomial rings\nFree modules over polynomial rings and submodules given by a finite generating set\nGroebner basis over a field\nFree/minimal resolutions\nSyzygy modules\nNemo.jl rings can be used as coefficient rings"
},

{
    "location": "index.html#Installation-1",
    "page": "Getting Started",
    "title": "Installation",
    "category": "section",
    "text": "To use Singular.jl we require Julia 0.6 or higher. Please see http://julialang.org/downloads for instructions on how to obtain julia for your system.At the Julia prompt simply typejulia> Pkg.clone(\"https://github.com/wbhart/Singular.jl\")\njulia> Pkg.build(\"Singular\")Note that Singular.jl depends on Cxx.jl which is not supported on every system."
},

{
    "location": "index.html#Quick-start-1",
    "page": "Getting Started",
    "title": "Quick start",
    "category": "section",
    "text": "Here is an example of using Singular.jljulia> using Singular\n\njulia> R, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n(Singular Polynomial Ring (QQ),(x,y),(dp(2),C), Singular.spoly{Singular.n_Q}[x, y])\n\njulia> I = Ideal(R, x^2 + 1, x*y + 1)\nSingular Ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x^2+1, x*y+1)\n\njulia> G = std(I)\nSingular Ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x-y, y^2+1)\n\njulia> Z = syz(G)\nSingular Module over Singular Polynomial Ring (QQ),(x,y),(dp(2),C), with Generators:\ny^2*gen(1)-x*gen(2)+y*gen(2)+gen(1)\n\njulia> F = fres(G, 0)\nSingular Resolution:\nR^1 <- R^2 <- R^1\n\njulia> F[1]\nSingular Module over Singular Polynomial Ring (QQ),(x,y),(dp(2),C), with Generators:\nx-y\ny^2+1"
},

{
    "location": "integer.html#",
    "page": "Integers",
    "title": "Integers",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "integer.html#Integers-1",
    "page": "Integers",
    "title": "Integers",
    "category": "section",
    "text": "The default integer type in Singular.jl is the Singular n_Z integer type.The associated ring of integers is represented by the constant parent object which can be constructed by a call to Singular.Integers().For convenience we defineZZ = Singular.Integers()so that integers can be constructed using ZZ. Note that this is the name of a specific parent object, not the name of its type.The types of the integer ring parent objects and elements of the associated rings of integers are given in the following table according to the library provding them.Library Element type Parent type\nSingular n_Z Singular.IntegersAll integer element types belong directly to the abstract type RingElem and all the integer ring parent object types belong to the abstract type Ring."
},

{
    "location": "integer.html#Integer-functionality-1",
    "page": "Integers",
    "title": "Integer functionality",
    "category": "section",
    "text": "Singular.jl integers implement the ring and possibly some parts of the Euclidean ring interfaces of AbstractAlgebra.jl.https://nemocas.github.io/AbstractAlgebra.jl/rings.htmlhttps://nemocas.github.io/AbstractAlgebra.jl/euclidean.htmlBelow, we describe the functionality that is specific to the Singular integer ring."
},

{
    "location": "integer.html#Constructors-1",
    "page": "Integers",
    "title": "Constructors",
    "category": "section",
    "text": "ZZ(n::Integer)Coerce a Julia integer value into the integer ring."
},

{
    "location": "integer.html#AbstractAlgebra.Generic.isunit-Tuple{Singular.n_Z}",
    "page": "Integers",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "Method",
    "text": "isunit(n::n_Z)\n\nReturn true if n is pm 1.\n\n\n\n"
},

{
    "location": "integer.html#Base.denominator-Tuple{Singular.n_Z}",
    "page": "Integers",
    "title": "Base.denominator",
    "category": "Method",
    "text": "denominator(n::n_Z)\n\nReturn the denominator of n (which will always be 1).\n\n\n\n"
},

{
    "location": "integer.html#Base.numerator-Tuple{Singular.n_Z}",
    "page": "Integers",
    "title": "Base.numerator",
    "category": "Method",
    "text": "numerator(n::n_Z)\n\nReturn the numerator of n (which is n itself).\n\n\n\n"
},

{
    "location": "integer.html#Base.abs-Tuple{Singular.n_Z}",
    "page": "Integers",
    "title": "Base.abs",
    "category": "Method",
    "text": "abs(n::n_Z)\n\nReturn the absolute value of n.\n\n\n\n"
},

{
    "location": "integer.html#Basic-manipulation-1",
    "page": "Integers",
    "title": "Basic manipulation",
    "category": "section",
    "text": "isunit(::n_Z)denominator(::n_Z)numerator(::n_Z)abs(::n_Z)Examplesa = ZZ(-12)\n\nisunit(a)\nn = numerator(a)\nd = denominator(a)\nc = abs(a)"
},

{
    "location": "integer.html#Euclidean-division-1",
    "page": "Integers",
    "title": "Euclidean division",
    "category": "section",
    "text": "Singular.jl provides a number of Euclidean division operations. Recall that for a dividend a and divisor b, we can write a = bq + r with 0 leq r  b. We call q the quotient and r the remainder.In the following table we list the division functions and their rounding behaviour. We also give the return value of the function, with q representing return of the quotient and r representing return of the remainder.Function Return Rounding\ndivrem(a::n_Z, b::n_Z) q, r towards zero\nrem(a::n_Z, b::n_Z) r towards zero\nmod(a::n_Z, b::n_Z) r downExamplesa = ZZ(-12)\nb = ZZ(5)\n\nq, r = divrem(a, b)\nr = mod(a, b)\nc = a % b"
},

{
    "location": "integer.html#Comparison-1",
    "page": "Integers",
    "title": "Comparison",
    "category": "section",
    "text": "Here is a list of the comparison functions implemented, with the understanding that isless provides all the usual comparison operators.Function\nisless(a::n_Z, b::n_Z)We also provide the following ad hoc comparisons which again provide all of the comparison operators mentioned above.Function\nisless(a::n_Z, b::Integer)\nisless(a::Integer, b::n_Z)Examplesa = ZZ(12)\nb = ZZ(3)\n\na < b\na != b\na > 4\n5 <= b"
},

{
    "location": "rational.html#",
    "page": "Rational field",
    "title": "Rational field",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "rational.html#Rational-field-1",
    "page": "Rational field",
    "title": "Rational field",
    "category": "section",
    "text": "Singular.jl provides rational numbers via Singular\'s n_Q type.There is a constant parent object representing the field of rationals, called QQ in Singular.jl. It is defined by QQ = Rationals(), which calls the constructor for the unique field of rationals in Singular."
},

{
    "location": "rational.html#Rational-functionality-1",
    "page": "Rational field",
    "title": "Rational functionality",
    "category": "section",
    "text": "The rationals in Singular.jl implement the Field interface defined by AbstractAlgebra.jl. They also implement the Fraction Field interface.https://nemocas.github.io/AbstractAlgebra.jl/fields.htmlhttps://nemocas.github.io/AbstractAlgebra.jl/fraction_fields.htmlWe describe here only the extra functionality provided by Singular that is not already described in those interfaces."
},

{
    "location": "rational.html#Constructors-1",
    "page": "Rational field",
    "title": "Constructors",
    "category": "section",
    "text": "In addition to the standard constructors required for the interfaces listed above, Singular.jl provides the following constructors.QQ(n::n_Z)\nQQ(n::fmpz)Construct a Singular rational from the given integer n."
},

{
    "location": "rational.html#Base.numerator-Tuple{Singular.n_Q}",
    "page": "Rational field",
    "title": "Base.numerator",
    "category": "Method",
    "text": "numerator(n::n_Q)\n\nReturn the numerator of the given fraction.\n\n\n\n"
},

{
    "location": "rational.html#Base.denominator-Tuple{Singular.n_Q}",
    "page": "Rational field",
    "title": "Base.denominator",
    "category": "Method",
    "text": "denominator(n::n_Q)\n\nReturn the denominator of the given fraction.\n\n\n\n"
},

{
    "location": "rational.html#Base.abs-Tuple{Singular.n_Q}",
    "page": "Rational field",
    "title": "Base.abs",
    "category": "Method",
    "text": "abs(n::n_Q)\n\nReturn the absolute value of the given fraction.\n\n\n\n"
},

{
    "location": "rational.html#Basic-manipulation-1",
    "page": "Rational field",
    "title": "Basic manipulation",
    "category": "section",
    "text": "numerator(::n_Q)denominator(::n_Q)abs(::n_Q)f = QQ(-12, 7)\n\nh = numerator(QQ)\nk = denominator(QQ)\nm = abs(f)"
},

{
    "location": "rational.html#Comparison-1",
    "page": "Rational field",
    "title": "Comparison",
    "category": "section",
    "text": "Here is a list of the comparison functions implemented, with the understanding that isless provides all the usual comparison operators.Function\nisless(a::n_Q, b::n_Q)We also provide the following ad hoc comparisons which again provide all of the comparison operators mentioned above.Function\nisless(a::n_Q, b::Integer)\nisless(a::Integer, b::n_Q)Examplesa = QQ(12, 7)\nb = QQ(-3, 5)\n\na > b\na != b\na > 1\n5 >= b"
},

{
    "location": "rational.html#Nemo.reconstruct-Tuple{Singular.n_Z,Singular.n_Z}",
    "page": "Rational field",
    "title": "Nemo.reconstruct",
    "category": "Method",
    "text": "reconstruct(x::n_Z, y::n_Z)\n\nGiven x modulo y, find rs such that x equiv rs pmody for values r and s satisfying the bound y  2(r + 1)(s + 1).\n\n\n\n"
},

{
    "location": "rational.html#Rational-reconstruction-1",
    "page": "Rational field",
    "title": "Rational reconstruction",
    "category": "section",
    "text": "reconstruct(::n_Z, ::n_Z)The following ad hoc versions of the same function also exist.reconstruct(::n_Z, ::Integer)\nreconstruct(::Integer, ::n_Z)Examplesq1 = reconstruct(ZZ(7), ZZ(3))\nq2 = reconstruct(ZZ(7), 5)"
},

{
    "location": "modn.html#",
    "page": "Integers mod n",
    "title": "Integers mod n",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "modn.html#Integers-mod-n-1",
    "page": "Integers mod n",
    "title": "Integers mod n",
    "category": "section",
    "text": "Integers mod n are implemented via the Singular n_Zn type for any positive modulus that can fit in a Julia Int.The associated ring of integers mod n is represented by a parent object which can be constructed by a call to the ResidueRing constructor.The types of the parent objects and elements of the associated rings of integers modulo n are given in the following table according to the library providing them.Library Element type Parent type\nSingular n_Zn Singular.N_ZnRingAll integer mod n element types belong directly to the abstract type RingElem and all the parent object types belong to the abstract type Ring."
},

{
    "location": "modn.html#Integer-mod-n-functionality-1",
    "page": "Integers mod n",
    "title": "Integer mod n functionality",
    "category": "section",
    "text": "Singular.jl integers modulo n implement the Ring and Residue Ring interfaces of AbstractAlgebra.jl.https://nemocas.github.io/AbstractAlgebra.jl/rings.htmlhttps://nemocas.github.io/AbstractAlgebra.jl/residue_rings.htmlParts of the Euclidean Ring interface may also be implemented, though Singular will report an error if division is meaningless (even after cancelling zero divisors).https://nemocas.github.io/AbstractAlgebra.jl/euclidean.htmlBelow, we describe the functionality that is specific to the Singular integers mod n ring and not already listed at the given links."
},

{
    "location": "modn.html#Constructors-1",
    "page": "Integers mod n",
    "title": "Constructors",
    "category": "section",
    "text": "Given a ring R of integers modulo n, we also have the following coercions in addition to the standard ones expected.R(n::n_Z)\nR(n::fmpz)Coerce a Singular or Flint integer value into the ring."
},

{
    "location": "modn.html#AbstractAlgebra.Generic.isunit-Tuple{Singular.n_Zn}",
    "page": "Integers mod n",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "Method",
    "text": "isunit(n::n_Zn)\n\nReturn true if the given value is a unit in the integers modulo n.\n\n\n\n"
},

{
    "location": "modn.html#AbstractAlgebra.Generic.characteristic-Tuple{Singular.N_ZnRing}",
    "page": "Integers mod n",
    "title": "AbstractAlgebra.Generic.characteristic",
    "category": "Method",
    "text": "characteristic(R::N_ZnRing)\n\nReturn the characteristic n of the ring.\n\n\n\n"
},

{
    "location": "modn.html#Basic-manipulation-1",
    "page": "Integers mod n",
    "title": "Basic manipulation",
    "category": "section",
    "text": "isunit(::n_Zn)Singular.characteristic(::N_ZnRing)ExamplesR = ResidueRing(ZZ, 26)\na = R(5)\n\nisunit(a)\nc = characteristic(R)"
},

{
    "location": "modp.html#",
    "page": "Integers mod p",
    "title": "Integers mod p",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "modp.html#Integers-mod-p-1",
    "page": "Integers mod p",
    "title": "Integers mod p",
    "category": "section",
    "text": "Integers mod a prime p are implemented via the Singular n_Zp type for any positive prime modulus less than 2^28.The associated field of integers mod p is represented by a parent object which can be constructed by a call to the Fp constructor.The types of the parent objects and elements of the associated fields of integers modulo p are given in the following table according to the library providing them.Library Element type Parent type\nSingular n_Zp Singular.N_ZpFieldAll integer mod p element types belong directly to the abstract type FieldElem and all the parent object types belong to the abstract type Field."
},

{
    "location": "modp.html#Integer-mod-p-functionality-1",
    "page": "Integers mod p",
    "title": "Integer mod p functionality",
    "category": "section",
    "text": "Singular.jl integers modulo p implement the Field and Residue Ring interfaces of AbstractAlgebra.jl.https://nemocas.github.io/AbstractAlgebra.jl/fields.htmlhttps://nemocas.github.io/AbstractAlgebra.jl/residue_rings.htmlBelow, we describe the functionality that is specific to the Singular integers mod p field and not already listed at the given links."
},

{
    "location": "modp.html#Constructors-1",
    "page": "Integers mod p",
    "title": "Constructors",
    "category": "section",
    "text": "The following constructors are available to create the field of integers modulo a prime p.Fp(p::Int; cached=true)Construct the field of integers modulo p. By default, the field is cached, so that all fields of integers modulo p have the same parent object. If this is not the desired behaviour, the cached paramater can be set to false. If p is not a prime or p is not in the range (0 2^28), an exception is raised.Given a field R of integers modulo p, we also have the following coercions in addition to the standard ones expected.R(n::n_Z)\nR(n::fmpz)Coerce a Singular or Flint integer value into the field."
},

{
    "location": "modp.html#AbstractAlgebra.Generic.isunit-Tuple{Singular.n_Zp}",
    "page": "Integers mod p",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "Method",
    "text": "isunit(n::n_Zp)\n\nReturn true if n is a unit in the field, i.e. nonzero.\n\n\n\n"
},

{
    "location": "modp.html#AbstractAlgebra.Generic.characteristic-Tuple{Singular.N_ZpField}",
    "page": "Integers mod p",
    "title": "AbstractAlgebra.Generic.characteristic",
    "category": "Method",
    "text": "characteristic(R::N_ZpField)\n\nReturn the characteristic of the field.\n\n\n\n"
},

{
    "location": "modp.html#Basic-manipulation-1",
    "page": "Integers mod p",
    "title": "Basic manipulation",
    "category": "section",
    "text": "isunit(::n_Zp)Singular.characteristic(::N_ZpField)ExamplesR = Fp(23)\na = R(5)\n\nisunit(a)\nc = characteristic(R)"
},

{
    "location": "modp.html#Conversions-1",
    "page": "Integers mod p",
    "title": "Conversions",
    "category": "section",
    "text": "Int(n::n_Zp)Lift the integer n modulo p to a Julia Int. The result is always in the range 0 p).ExamplesR = Fp(23)\na = R(5)\n\nb = Int(a)"
},

{
    "location": "GF.html#",
    "page": "Finite fields",
    "title": "Finite fields",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "GF.html#Finite-fields-1",
    "page": "Finite fields",
    "title": "Finite fields",
    "category": "section",
    "text": "Finite fields are implemented via the Singular n_GF type for any characteristic and degree contained in the Singular Conway tables.The associated finite field is represented by a parent object which can be constructed by a call to the FiniteField constructor.The types of the parent objects and elements of the associated finite fields are given in the following table according to the library providing them.Library Element type Parent type\nSingular n_GF Singular.N_GFieldAll finite field element types belong directly to the abstract type FieldElem and all the parent object types belong to the abstract type Field."
},

{
    "location": "GF.html#Finite-field-functionality-1",
    "page": "Finite fields",
    "title": "Finite field functionality",
    "category": "section",
    "text": "Singular.jl finite fields implement the Field interface of AbstractAlgebra.jl.https://nemocas.github.io/AbstractAlgebra.jl/fields.htmlBelow, we describe the functionality that is specific to Singular finite field and not already listed at the given link."
},

{
    "location": "GF.html#Singular.FiniteField-Tuple{Int64,Int64,String}",
    "page": "Finite fields",
    "title": "Singular.FiniteField",
    "category": "Method",
    "text": "FiniteField(p::Int, n::Int, S::String; cached=true)\n\nReturns a tuple K, a consisting of a finite field K of characteristic p and degree n, and its generator a. The string used to print the generator is given by S. If the finite field is not listed in the Conway tables included in Singular, an error will be raised. By default, finite fields are cached globally, so that there is only one finite field in the system with given characteristic, degree and string. If this is not the desired behaviour, one can pass false for the optional cached parameter.\n\n\n\n"
},

{
    "location": "GF.html#Constructors-1",
    "page": "Finite fields",
    "title": "Constructors",
    "category": "section",
    "text": "The following constructors are available to create finite fields and their elements.Singular.FiniteField(::Int, ::Int, ::String; ::Bool)Given a finite field R, we also have the following coercions in addition to the standard ones expected.R(n::fmpz)Coerce a Flint integer value into the field."
},

{
    "location": "GF.html#AbstractAlgebra.Generic.degree-Tuple{Singular.N_GField}",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.degree",
    "category": "Method",
    "text": "degree(R::N_GField)\n\nReturn the degree of the field as an extension of mathbbF_p.\n\n\n\n"
},

{
    "location": "GF.html#AbstractAlgebra.Generic.characteristic-Tuple{Singular.N_GField}",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.characteristic",
    "category": "Method",
    "text": "characteristic(R::N_GField)\n\nReturn the characteristic of the field.\n\n\n\n"
},

{
    "location": "GF.html#AbstractAlgebra.Generic.isunit-Tuple{Singular.n_GF}",
    "page": "Finite fields",
    "title": "AbstractAlgebra.Generic.isunit",
    "category": "Method",
    "text": "isunit(n::n_GF)\n\nReturn true if n is a unit in the field, i.e. nonzero.\n\n\n\n"
},

{
    "location": "GF.html#Basic-manipulation-1",
    "page": "Finite fields",
    "title": "Basic manipulation",
    "category": "section",
    "text": "Singular.degree(::N_GField)Singular.characteristic(::N_GField)isunit(::n_GF)ExamplesR = FiniteField(7, 2)\na = R(5)\n\nisunit(a)\nc = characteristic(R)\nd = degree(R)"
},

{
    "location": "nemo.html#",
    "page": "Nemo rings and fields",
    "title": "Nemo rings and fields",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "nemo.html#Nemo-rings-and-fields-1",
    "page": "Nemo rings and fields",
    "title": "Nemo rings and fields",
    "category": "section",
    "text": "Any type that satisfies AbstractAlgebra.jl Ring or Field interface, such as all Nemo ring and field types, can be used as coefficient rings in Singular.jl. Theses are implemented via the Singular n_unknown type.The associated ring/field is represented by a parent object which can be constructed by a call to the CoefficientRing constructor. In practice, however, this constructor is only used internally, and Nemo rings and fields work directly as Singular coefficient rings, and all the coercions and ad hoc functions that one would expect to be present are implemented.All of the Singular polynomial arithmetic should work for any Nemo ring and everything, including ideals, modules, standard basis, syzygies, resolutions, etc., should work with any Nemo field.The types of the parent objects and elements of the associated foreign rings are given in the following table according to the library providing them.Library Element type Parent type\nSingular n_unknown{T} Singular.CoefficientRing{T}These types are parameterised with the element type of the given Nemo/AbstractAlgebra element type.The Singular.jl n_unknown types belong directly to the abstract type RingElem and their parent object types belong to the abstract type Ring.Specialised efficient wrappers exist for certain Nemo coefficient ring types."
},

{
    "location": "nemo.html#Nemo-ring-functionality-1",
    "page": "Nemo rings and fields",
    "title": "Nemo ring functionality",
    "category": "section",
    "text": "Singular.jl foreign ring types implement the Ring interface and possibly the Field interface of AbstractAlgebra.jl.https://nemocas.github.io/AbstractAlgebra.jl/rings.htmlhttps://nemocas.github.io/AbstractAlgebra.jl/fields.htmlParts of the Euclidean Ring interface may also be implemented, though Singular will report an error if division is meaningless (even after cancelling zero divisors).https://nemocas.github.io/AbstractAlgebra.jl/euclidean.htmlBelow, we describe the functionality that is specific to the Singular foreign ring interface that is not already listed at the given links."
},

{
    "location": "nemo.html#Constructors-1",
    "page": "Nemo rings and fields",
    "title": "Constructors",
    "category": "section",
    "text": "Given an AbstractAlgebra compatible ring R, e.g. a Nemo ring, we have the following constructor, which returns the associated Singular.jl coefficient ring..CoefficientRing(R::Ring)If there are generators to be coerced from Nemo/AbstractAlgebra into corresponding elements, the Singular.jl coefficient ring can be used to coerce them to a Singular n_unknown element.ExamplesR, x = Nemo.PolynomialRing(ZZ, \"x\")\nS = CoefficientRing(R)\nt = S(x)Note that it is unlikely that a user directly needs to construct the Singular coefficient ring from a Nemo ring, since the Singular.jl constructors are designed to accept Nemo coefficient rings directly. Singular.jl automatically constructs the required Singular coefficient ring and makes use of it."
},

{
    "location": "polynomial.html#",
    "page": "Multivariate polynomials",
    "title": "Multivariate polynomials",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "polynomial.html#Multivariate-polynomials-1",
    "page": "Multivariate polynomials",
    "title": "Multivariate polynomials",
    "category": "section",
    "text": "Singular.jl allows the creation of multivariate polynomials over any of the coefficient rings described above.The default multivariate polynomial type in Singular.jl is the Singular spoly type.The associated polynomial ring is represented by a parent object which can be constructed by a call to the PolynomialRing constructor.The types of the polynomial ring parent objects and elements thereof are given in the following table according to the library provding them.Library Element type Parent type\nSingular spoly{T} Singular.PolyRing{T}These types are parameterised by the type of elements in the coefficient ring of the polynomials.All polynomial types belong directly to the abstract type RingElem and all the polynomial ring parent object types belong to the abstract type Ring."
},

{
    "location": "polynomial.html#Multivariate-polynomial-functionality-1",
    "page": "Multivariate polynomials",
    "title": "Multivariate polynomial functionality",
    "category": "section",
    "text": "Singular.jl polynomials implement the Multivariate Polynomial Ring interface of AbstractAlgebra.jl.https://nemocas.github.io/AbstractAlgebra.jl/mpolynomial_rings.htmlIn particular, Singular polynomials are sparse distributed, but do not have random access. Instead, they implement iterator access to terms. This is due to their storage in a linked list, for efficient implementation of Groebner basis algorithms.Some polynomial rings may also implement part of the Euclidean Ring interface, where this is appropriate.https://nemocas.github.io/AbstractAlgebra.jl/euclidean.htmlBelow, we describe the functionality that is specific to the Singular multivariate  polynomials that is not documented in the general multivariate interface."
},

{
    "location": "polynomial.html#Constructors-1",
    "page": "Multivariate polynomials",
    "title": "Constructors",
    "category": "section",
    "text": "PolynomialRing(R::Union{Ring, Field}, s::Array{String, 1};\n   cached::Bool = true, ordering::Symbol = :degrevlex,\n      ordering2::Symbol = :comp1min, degree_bound::Int = 0)Returns a tuple, S x consisting of a multivariate polynomial ring S and an array of variables (from which polynomials can be constructed). The ring R must be a valid Singular coefficient ring, or any Nemo/AbstractAlgebra coefficient ring. The array s must be a list of strings corresponding to how the variables will be printed. By default, there will only be one Singular polynomial ring in the system for each combination of coefficient ring, list of variable names, ordering and degree bound. This is accomplished by making use of a global cache. If this is not the desired behaviour, false can be passed to the optional argument cached. Two orderings can be specified, one for term ordering of the polynomials, and another for ordering of module components. They can occur in either order, the first taking precedence over the other, when the polynomials are used to represent module generators. If either is not specified, the indicated default is used.The options for polynomial term ordering are the symbols, :lex, :deglex and :degrevlex, and the options for module component ordering are comp1min and comp1max.If one has an a priori bound on the degree in each variable of a polynomial (including for all intermediate computations in this ring), one can specify it using the degree_bound optional parameter. Singular may then be able to use a more efficient representation internally, which will use less memory and allow for faster arithmetic. By default, Singular uses a bound of 16 bits internally for the exponent of each variable, however this is a signed value, so that the default is for nonnegative exponents that fit in 15 bits.Note that internally, Singular may use a higher bound than specified, if it will not increase the amount of storage required.ExamplesR, (x, y, z) = PolynomialRing(ZZ, [\"x\", \"y\", \"z\"])\n\nS, vars = PolynomialRing(QQ, [\"x\", \"y\"]; ordering=:deglex)\n\nT, x = PolynomialRing(ZZ, [\"x$i\" for i in 1:5];\n       ordering=:comp1max, ordering2=:degrevlex, degree_bound=5)See also the convenience macros below for simple use cases."
},

{
    "location": "polynomial.html#Polynomial-ring-macros-1",
    "page": "Multivariate polynomials",
    "title": "Polynomial ring macros",
    "category": "section",
    "text": "For convenience, we provide some macros for constructing polynomial rings and injecting the variables into scope. These are easier to use, but have some limitations. In particular, they can only be used at the top level by the user, and cannot be used programmatically or in library code (it is not possible to inject an arbitrary number of variables into scope inside a function).The macros are designed for simple use cases, and do not offer the full power of the most general constructor above.@PolynomialRing(R, s, n, o)Given a coefficient ring R, a root variable name, e.g. \"x\", a number of variable n and a polynomial term ordering o, create the variables x1, x2, ..., xn and inject them into scope, and return the corresponding polynomial ring S.@PolynomialRing(R, s, n)As per the previous macro, with a default of :degrevlex for the polynomial term ordering.ExamplesS = @PolynomialRing(ZZ, \"x\", 5, :deglex)\n\nT = @PolynomialRing(QQ, \"y\", 10)"
},

{
    "location": "polynomial.html#Singular.ngens-Tuple{Singular.PolyRing}",
    "page": "Multivariate polynomials",
    "title": "Singular.ngens",
    "category": "Method",
    "text": "ngens(R::PolyRing)\n\nReturn the number of variables in the given polynomial ring.\n\n\n\n"
},

{
    "location": "polynomial.html#Singular.has_global_ordering-Tuple{Singular.PolyRing}",
    "page": "Multivariate polynomials",
    "title": "Singular.has_global_ordering",
    "category": "Method",
    "text": "has_global_ordering(R::PolyRing)\n\nReturn true if the given polynomial has a global ordering, i.e. if 1  x for each variable x in the ring. This include :lex, :deglex and :degrevlex orderings..\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.characteristic-Tuple{Singular.PolyRing}",
    "page": "Multivariate polynomials",
    "title": "AbstractAlgebra.Generic.characteristic",
    "category": "Method",
    "text": "characteristic(R::PolyRing)\n\nReturn the characteristic of the polynomial ring, i.e. the characteristic of the coefficient ring.\n\n\n\n"
},

{
    "location": "polynomial.html#Singular.degree_bound-Tuple{Singular.PolyRing}",
    "page": "Multivariate polynomials",
    "title": "Singular.degree_bound",
    "category": "Method",
    "text": "degree_bound(R::PolyRing)\n\nReturn the internal degree bound in each variable, enforced by Singular. This is the largest positive value any degree can have before an overflow will occur. This internal bound may be higher than the bound requested by the user via the  degree_bound parameter of the PolynomialRing constructor.\n\n\n\n"
},

{
    "location": "polynomial.html#Singular.lead_exponent-Tuple{Singular.spoly}",
    "page": "Multivariate polynomials",
    "title": "Singular.lead_exponent",
    "category": "Method",
    "text": "lead_exponent(p::spoly)\n\nReturn the exponent vector of the leading term of the given polynomial. The return value is a Julia 1-dimensional array giving the exponent for each variable of the leading term.\n\n\n\n"
},

{
    "location": "polynomial.html#Basic-manipulation-1",
    "page": "Multivariate polynomials",
    "title": "Basic manipulation",
    "category": "section",
    "text": "ngens(::PolyRing)has_global_ordering(R::PolyRing)Singular.characteristic(R::PolyRing)degree_bound(R::PolyRing)lead_exponent(p::spoly)ExamplesR = @PolynomialRing(ZZ, \"x\", 3)\n\nn = ngens(R)\nhas_global_ordering(R) == true\nc = characteristic(R)\nL = degree_bound(R)\nexps = lead_exponent(x1*x2 + 3x1*x2^2 + x3 + 2)"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.primpart-Tuple{Singular.spoly}",
    "page": "Multivariate polynomials",
    "title": "AbstractAlgebra.Generic.primpart",
    "category": "Method",
    "text": "primpart(x::spoly)\n\nReturn the primitive part of the polynomial, i.e. the polynomial divided by the GCD of its coefficients.\n\n\n\n"
},

{
    "location": "polynomial.html#AbstractAlgebra.Generic.content-Tuple{Singular.spoly}",
    "page": "Multivariate polynomials",
    "title": "AbstractAlgebra.Generic.content",
    "category": "Method",
    "text": "content(x::spoly)\n\nReturn the content of the polynomial, i.e. the GCD of its coefficients.\n\n\n\n"
},

{
    "location": "polynomial.html#Content-and-primitive-part-1",
    "page": "Multivariate polynomials",
    "title": "Content and primitive part",
    "category": "section",
    "text": "When coefficient rings have a meaningful GCD function, the following functions are available.Singular.primpart(x::spoly)Singular.content(x::spoly)ExamplesR = @PolynomialRing(ZZ, \"x\", 2)\n\nf = 3x1^2 + 3x1*x2 + 6x2^2\n\np = primpart(f)\nc = content(f)"
},

{
    "location": "ideal.html#",
    "page": "Ideals",
    "title": "Ideals",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "ideal.html#Ideals-1",
    "page": "Ideals",
    "title": "Ideals",
    "category": "section",
    "text": "Singular.jl allows the creation of ideals over a Singular polynomial ring. These are internally stored as a list of (polynomial) generators. This list of generators can also  have the property of being a Groebner basis.The default ideal type in Singular.jl is the Singular sideal type.Ideals objects have a parent object which represents the set of ideals they belong to, the data for which is given by the polynomial ring their generators belong to.The types of ideals and associated parent objects are given in the following table according to the library provding them.Library Element type Parent type\nSingular sideal{T} Singular.IdealSet{T}These types are parameterised by the type of elements of the polynomial ring over which the ideals are defined.All ideal types belong directly to the abstract type Module{T} and all the ideal set parent object types belong to the abstract type Set."
},

{
    "location": "ideal.html#Ideal-functionality-1",
    "page": "Ideals",
    "title": "Ideal functionality",
    "category": "section",
    "text": "Singular.jl ideals implement standard operations one would expect on modules. These include:Operations common to all AbstractAlgebra objects, such as parent, base_ring, elem_type, parent_type, parent, deepcopy, etc.\nAdditionAlso implements is the following operations one expects for ideals:Multiplication\nPoweringBelow, we describe all of the functionality for Singular.jl ideals that is not included in this list of basic operations."
},

{
    "location": "ideal.html#Constructors-1",
    "page": "Ideals",
    "title": "Constructors",
    "category": "section",
    "text": "Given a Singular polynomial ring R, the following constructors are available for creating ideals.Ideal(R::PolyRing{T}, ids::spoly{T}...) where T <: Nemo.RingElem\nIdeal(R::PolyRing{T}, ids::Array{spoly{T}, 1}) where T <: Nemo.RingElemConstruct the ideal over the polynomial ring R whose (polynomial) generators are given  by the given parameter list or array of polynomials, respectively. The list may be empty, resulting in the zero ideal.ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nI1 = Ideal(R, x*y + 1, x^2)\nI2 = Ideal(R, [x*y + 1, x^2])"
},

{
    "location": "ideal.html#Singular.ngens-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "Singular.ngens",
    "category": "Method",
    "text": "ngens(I::sideal)\n\nReturn the number of generators in the internal representation of the ideal I.\n\n\n\n"
},

{
    "location": "ideal.html#Base.iszero-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "Base.iszero",
    "category": "Method",
    "text": "iszero(I::sideal)\n\nReturn true if the given ideal is algebraically the zero ideal.\n\n\n\n"
},

{
    "location": "ideal.html#Singular.iszerodim-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "Singular.iszerodim",
    "category": "Method",
    "text": "iszerodim(I::sideal)\n\nReturn true if the given ideal is zero dimensional, i.e. the Krull dimension of RI is zero, where R is the polynomial ring over which I is an ideal..\n\n\n\n"
},

{
    "location": "ideal.html#AbstractAlgebra.Generic.isconstant-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "AbstractAlgebra.Generic.isconstant",
    "category": "Method",
    "text": "isconstant(I::sideal)\n\nReturn true if the given ideal is a constant ideal, i.e. generated by constants in the polynomial ring over which it is an ideal.\n\n\n\n"
},

{
    "location": "ideal.html#Singular.isvar_generated-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "Singular.isvar_generated",
    "category": "Method",
    "text": "isvar_generated(I::sideal)\n\nReturn true if each generator in the representation of the ideal I is a generator of the polynomial ring, i.e. a variable.\n\n\n\n"
},

{
    "location": "ideal.html#Base.LinAlg.normalize!-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "Base.LinAlg.normalize!",
    "category": "Method",
    "text": "normalize!(I::sideal)\n\nNormalize the polynomial generators of the ideal I in-place. This means to reduce their coefficients to lowest terms. In most cases this does nothing, but if the coefficient ring were the rational numbers for example, the coefficients of the polynomials would be reduced to lowest terms.\n\n\n\n"
},

{
    "location": "ideal.html#Basic-manipulation-1",
    "page": "Ideals",
    "title": "Basic manipulation",
    "category": "section",
    "text": "ngens(::sideal)Singular.jl overloads the setindex! and getindex functions so that one can access the generators of an ideal using array notation.I[n::Int]iszero(::sideal)iszerodim(::sideal)isconstant(::sideal)isvar_generated(::sideal)normalize!(::sideal)ExamplesR, (x, y) = PolynomialRing(ZZ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2 + 1, x*y)\n\nn = ngens(I)\np = I[1]\nI[1] = 2x + y^2\nisconstant(I) == false\nisvar_generated(I) == false\niszerodim(I) == false"
},

{
    "location": "ideal.html#Base.contains-Union{Tuple{Singular.sideal{T},Singular.sideal{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Base.contains",
    "category": "Method",
    "text": "contains{T <: AbstractAlgebra.RingElem}(I::sideal{T}, J::sideal{T})\n\nReturns true if the ideal I contains the ideal J. This will be expensive if I is not a Groebner ideal, since its standard basis must be computed.\n\n\n\n"
},

{
    "location": "ideal.html#Containment-1",
    "page": "Ideals",
    "title": "Containment",
    "category": "section",
    "text": "contains{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})ExamplesR, (x , y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2 + 1, x*y)\nJ = Ideal(R, x^2 + 1)\n\ncontains(I, J) == true"
},

{
    "location": "ideal.html#Base.isequal-Union{Tuple{Singular.sideal{T},Singular.sideal{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Base.isequal",
    "category": "Method",
    "text": "isequal{T <: AbstractAlgebra.RingElem}(I1::sideal{T}, I2::sideal{T})\n\nReturn true if the given ideals have the same generators in the same order. Note that two algebraically equal ideals with different generators will return false.\n\n\n\n"
},

{
    "location": "ideal.html#Singular.equal-Union{Tuple{Singular.sideal{T},Singular.sideal{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Singular.equal",
    "category": "Method",
    "text": "equal(I1::sideal{T}, I2::sideal{T}) where T <: AbstractAlgebra.RingElem\n\nReturn true if the two ideals are contained in each other, i.e. are the same ideal mathematically. This function should be called only as a last resort; it is exceptionally expensive to test equality of ideals! Do not define == as an alias for this function!\n\n\n\n"
},

{
    "location": "ideal.html#Comparison-1",
    "page": "Ideals",
    "title": "Comparison",
    "category": "section",
    "text": "Checking whether two ideals are algebraically equal is very expensive, as it usually requires computing Groebner bases. Therefore we do not overload the == operator for ideals. Instead we have the following two functions.isequal{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})equal{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})ExamplesR, (x , y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2 + 1, x*y)\nJ = Ideal(R, x^2 + x*y + 1, x^2 - x*y + 1)\n\nisequal(I, J) == false\nequal(I, J) == true"
},

{
    "location": "ideal.html#Singular.intersection-Union{Tuple{Singular.sideal{T},Singular.sideal{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Singular.intersection",
    "category": "Method",
    "text": "intersection{T <: Nemo.RingElem}(I::sideal{T}, J::sideal{T})\n\nReturns the intersection of the two given ideals.\n\n\n\n"
},

{
    "location": "ideal.html#Intersection-1",
    "page": "Ideals",
    "title": "Intersection",
    "category": "section",
    "text": "intersection{T <: Nemo.RingElem}(::sideal{T}, ::sideal{T})ExamplesR, (x , y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2 + 1, x*y)\nJ = Ideal(R, x^2 + x*y + 1, x^2 - x*y + 1)\n\nV = intersection(I, J)"
},

{
    "location": "ideal.html#Singular.quotient-Union{Tuple{Singular.sideal{T},Singular.sideal{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Singular.quotient",
    "category": "Method",
    "text": "quotient{T <: Nemo.RingElem}(I::sideal{T}, J::sideal{T})\n\nReturns the quotient of the two given ideals. Recall that the ideal quotient (IJ) over a polynomial ring R is defined by r in R  rJ subseteq I. \n\n\n\n"
},

{
    "location": "ideal.html#Quotient-1",
    "page": "Ideals",
    "title": "Quotient",
    "category": "section",
    "text": "quotient{T <: Nemo.RingElem}(::sideal{T}, ::sideal{T})ExamplesR, (x , y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2 + 1, x*y)\nJ = Ideal(R, x + y)\n\nV = quotient(I, J)"
},

{
    "location": "ideal.html#AbstractAlgebra.Generic.lead-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "AbstractAlgebra.Generic.lead",
    "category": "Method",
    "text": "lead(I::sideal)\n\nReturn the ideal generated by the leading terms of the polynomials generating I.\n\n\n\n"
},

{
    "location": "ideal.html#Leading-terms-1",
    "page": "Ideals",
    "title": "Leading terms",
    "category": "section",
    "text": "lead(::sideal)ExamplesR, (x , y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2 + 1, x*y)\n\nV = lead(I)"
},

{
    "location": "ideal.html#Singular.saturation-Union{Tuple{Singular.sideal{T},Singular.sideal{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Singular.saturation",
    "category": "Method",
    "text": "saturation{T <: Nemo.RingElem}(I::sideal{T}, J::sideal{T})\n\nReturns the saturation of the ideal I with respect to J, i.e. returns the quotient ideal (IJ^infty).\n\n\n\n"
},

{
    "location": "ideal.html#Saturation-1",
    "page": "Ideals",
    "title": "Saturation",
    "category": "section",
    "text": "saturation{T <: Nemo.RingElem}(::sideal{T}, ::sideal{T})ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, (x^2 + x*y + 1)*(2y^2+1)^3, (2y^2 + 3)*(2y^2+1)^2)\nJ = Ideal(R, 2y^2 + 1)\n\nS = saturation(I, J)"
},

{
    "location": "ideal.html#Base.std-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "Base.std",
    "category": "Method",
    "text": "std(I::sideal; complete_reduction::Bool=false)\n\nCompute a Groebner basis for the ideal I. Note that without complete_reduction set to true, the generators of the Groebner basis only have unique leading terms (up to permutation and multiplication by constants). If complete_reduction is set to true (and the ordering is a global ordering) then the Groebner basis is unique.\n\n\n\n"
},

{
    "location": "ideal.html#Singular.satstd-Union{Tuple{Singular.sideal{T},Singular.sideal{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Singular.satstd",
    "category": "Method",
    "text": "satstd{T <: AbstractAlgebra.RingElem}(I::sideal{T}, J::sideal{T})\n\nGiven an ideal J generated by variables, computes a standard basis of saturation(I, J). This is accomplished by dividing polynomials that occur throughout the std computation by variables occuring in J, where possible. Thus the result can be obtained faster than by first computing the saturation and then the standard basis.\n\n\n\n"
},

{
    "location": "ideal.html#Standard-basis-1",
    "page": "Ideals",
    "title": "Standard basis",
    "category": "section",
    "text": "std(::sideal; ::Bool)satstd{T <: AbstractAlgebra.RingElem}(::sideal{T}, ::sideal{T})ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2 + x*y + 1, 2y^2 + 3)\nJ = Ideal(R, 2*y^2 + 3, x^2 + x*y + 1)\n\nA = std(I)\n\nR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, (x*y + 1)*(2x^2*y^2 + x*y - 2) + 2x*y^2 + x, 2x*y + 1)\nJ = Ideal(R, x)\n\nB = satstd(I, J)"
},

{
    "location": "ideal.html#Base.reduce-Tuple{Singular.sideal,Singular.sideal}",
    "page": "Ideals",
    "title": "Base.reduce",
    "category": "Method",
    "text": "reduce(I::sideal, G::sideal)\n\nReturn an ideal whose generators are the generators of I reduced by the ideal G. The ideal G is required to be a Groebner basis. The returned ideal will have the same number of generators as I, even if they are zero.\n\n\n\n"
},

{
    "location": "ideal.html#Base.reduce-Tuple{Singular.spoly,Singular.sideal}",
    "page": "Ideals",
    "title": "Base.reduce",
    "category": "Method",
    "text": "reduce(p::spoly, G::sideal)\n\nReturn the polynomial which is p reduced by the polynomials generating G. It is assumed that G is a Groebner basis.\n\n\n\n"
},

{
    "location": "ideal.html#Reduction-1",
    "page": "Ideals",
    "title": "Reduction",
    "category": "section",
    "text": "reduce(::sideal, ::sideal)reduce(::spoly, ::sideal)ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nf = x^2*y + 2y + 1\ng = y^2 + 1\n\nI = Ideal(R, (x^2 + 1)*f + (x + y)*g + x + 1, (2y^2 + x)*f + y)\nJ = std(Ideal(R, f, g))\n\nV = reduce(I, J)\n\nh1 = (x^2 + 1)*f + (x + y)*g + x + 1\n\nh2 = reduce(h1, J)"
},

{
    "location": "ideal.html#Singular.eliminate-Tuple{Singular.sideal,Vararg{Singular.spoly,N} where N}",
    "page": "Ideals",
    "title": "Singular.eliminate",
    "category": "Method",
    "text": "eliminate(I::sideal, polys::spoly...)\n\nGiven a list of polynomials which are variables, construct the ideal corresponding geometrically to the projection of the variety given by the ideal I where those variables have been eliminated.\n\n\n\n"
},

{
    "location": "ideal.html#Elimination-1",
    "page": "Ideals",
    "title": "Elimination",
    "category": "section",
    "text": "eliminate(::sideal, ::spoly...)ExamplesR, (x, y, t) = PolynomialRing(QQ, [\"x\", \"y\", \"t\"])\n\nI = Ideal(R, x - t^2, y - t^3)\n\nJ = eliminate(I, t)"
},

{
    "location": "ideal.html#Singular.syz-Tuple{Singular.sideal}",
    "page": "Ideals",
    "title": "Singular.syz",
    "category": "Method",
    "text": "syz(I::sideal)\n\nCompute the module of syzygies of the ideal.\n\n\n\n"
},

{
    "location": "ideal.html#Syzygies-1",
    "page": "Ideals",
    "title": "Syzygies",
    "category": "section",
    "text": "syz(::sideal)ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2*y + 2y + 1, y^2 + 1)\n\nF = syz(I)\n\nM = Singular.Matrix(I)\nN = Singular.Matrix(F)\n\n# check they are actually syzygies\niszero(M*N)"
},

{
    "location": "ideal.html#Singular.fres-Union{Tuple{Singular.sideal{T},Int64,String}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Singular.fres",
    "category": "Method",
    "text": " fres{T <: Nemo.RingElem}(id::sideal{T}, max_length::Int,\n  method::String=\"complete\")\n\nCompute a free resolution of the given ideal up to the maximum given length. The ideal must be over a polynomial ring over a field, and a Groebner basis. The possible methods are \"complete\", \"frame\", \"extended frame\" and \"single module\". The result is given as a resolution, whose i-th entry is the syzygy module of the previous module, starting with the given ideal. The max_length can be set to 0 if the full free resolution is required.\n\n\n\n"
},

{
    "location": "ideal.html#Singular.sres-Union{Tuple{Singular.sideal{T},Int64}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Ideals",
    "title": "Singular.sres",
    "category": "Method",
    "text": " sres{T <: Nemo.RingElem}(id::sideal{T}, max_length::Int)\n\nCompute a (free) Schreyer resolution of the given ideal up to the maximum given length. The ideal must be over a polynomial ring over a field, and a Groebner basis. The result is given as a resolution, whose i-th entry is the syzygy module of the previous module, starting with the given ideal. The max_length can be set to 0 if the full free resolution is required.\n\n\n\n"
},

{
    "location": "ideal.html#Free-resolutions-1",
    "page": "Ideals",
    "title": "Free resolutions",
    "category": "section",
    "text": "fres{T <: Nemo.RingElem}(::sideal{T}, ::Int, ::String)sres{T <: Nemo.RingElem}(::sideal{T}, ::Int)ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x^2*y + 2y + 1, y^2 + 1)\n\nF1 = fres(std(I), 0)\nF2 = sres(std(I), 2)"
},

{
    "location": "module.html#",
    "page": "Finitely generated modules",
    "title": "Finitely generated modules",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "module.html#Finitely-generated-modules-1",
    "page": "Finitely generated modules",
    "title": "Finitely generated modules",
    "category": "section",
    "text": "Singular.jl allows the creation of submodules of a free module over a Singular polynomial ring, given by a finite generating set. These are internally stored as a list of elements of a free module over a polynomial ring R. This list of generators can also  have the property of being a Groebner basis.The default finitely generated module type in Singular.jl is the Singular smodule type.Module objects have a parent object which represents the class of R-modules they belong to, the data for which is given by the polynomial ring R over which the modules are defined.The types of modules and associated parent objects are given in the following table according to the library provding them.Library Element type Parent type\nSingular smodule{T} Singular.ModuleClass{T}These types are parameterised by the type of elements in the polynomial ring R.All module types belong directly to the abstract type Module{T} and all the module class parent object types belong to the abstract type Set."
},

{
    "location": "module.html#Module-functionality-1",
    "page": "Finitely generated modules",
    "title": "Module functionality",
    "category": "section",
    "text": "Singular.jl modules implement standard operations one would expect on modules. These include:Operations common to all AbstractAlgebra objects, such as parent, base_ring, elem_type, parent_type, parent, deepcopy, etc.Below, we describe all of the functionality for Singular.jl modules that is not included in this list of basic operations."
},

{
    "location": "module.html#Constructors-1",
    "page": "Finitely generated modules",
    "title": "Constructors",
    "category": "section",
    "text": "Given a Singular polynomial ring R, the following constructors are available for creating modules.Module{T <: Nemo.RingElem}(R::PolyRing{T}, vecs::svector{spoly{T}}...)Construct the module over the polynomial ring R whose generators are given  by the given parameter list of vectors (of length n), each component of which is a polynomial. These vectors represent elements of the free module R^n.Note that Module must be prepended with the package name Singular to disambiguate from Base.Module.ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nv1 = vector(R, x + 1, x*y + 1, y)\nv2 = vector(R, x^2 + 1, 2x + 3y, x)\n\nM = Singular.Module(R, v1, v2)"
},

{
    "location": "module.html#Singular.ngens-Tuple{Singular.smodule}",
    "page": "Finitely generated modules",
    "title": "Singular.ngens",
    "category": "Method",
    "text": "ngens(I::smodule)\n\nReturn the number of generators in the current representation of the module (as a list of vectors).\n\n\n\n"
},

{
    "location": "module.html#Base.LinAlg.rank-Tuple{Singular.smodule}",
    "page": "Finitely generated modules",
    "title": "Base.LinAlg.rank",
    "category": "Method",
    "text": "rank(I::smodule)\n\nReturn the rank n of the ambient space R^n of which this module is a submodule.\n\n\n\n"
},

{
    "location": "module.html#Base.iszero-Tuple{Singular.smodule}",
    "page": "Finitely generated modules",
    "title": "Base.iszero",
    "category": "Method",
    "text": "iszero(p::smodule)\n\nReturn true if this is algebraically the zero module.\n\n\n\n"
},

{
    "location": "module.html#Basic-manipulation-1",
    "page": "Finitely generated modules",
    "title": "Basic manipulation",
    "category": "section",
    "text": "ngens(::smodule)rank(::smodule)Singular.jl overloads the setindex! and getindex functions so that one can access the generators of a module using array notation. Each entry is a vector in R^n.M[n::Int]iszero(::smodule)ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nv1 = vector(R, x + 1, x*y + 1, y)\nv2 = vector(R, x^2 + 1, 2x + 3y, x)\n\nM = Singular.Module(R, v1, v2)\n\niszero(M) == false\nM[1] == v1\nn = rank(M)\nd = ngens(M)"
},

{
    "location": "module.html#Base.std-Tuple{Singular.smodule}",
    "page": "Finitely generated modules",
    "title": "Base.std",
    "category": "Method",
    "text": "std(I::smodule; complete_reduction::Bool=false)\n\nCompute the Groebner basis of the module I. If complete_reduction is set to true, the result is unique, up to permutation of the generators and multiplication by constants. If not, only the leading terms are unique (up to permutation of the generators and multiplication by constants, of course). Presently the polynomial ring used must be over a field or over the Singular integers.\n\n\n\n"
},

{
    "location": "module.html#Standard-basis-1",
    "page": "Finitely generated modules",
    "title": "Standard basis",
    "category": "section",
    "text": "std(::smodule; ::Bool)ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nv1 = vector(R, x + 1, x*y + 1, y)\nv2 = vector(R, x^2 + 1, 2x + 3y, x)\nv3 = x*v1 + y*v2 + vector(R, x, y + 1, y^2)\n\nM = Singular.Module(R, v1, v2, v3)\n\nG = std(M; complete_reduction=true)"
},

{
    "location": "module.html#Singular.syz-Tuple{Singular.smodule}",
    "page": "Finitely generated modules",
    "title": "Singular.syz",
    "category": "Method",
    "text": "syz(M::smodule)\n\nCompute the module of syzygies of the given module. This will be given as a set of generators in an ambient space R^n, where n is the number of generators in M.\n\n\n\n"
},

{
    "location": "module.html#Syzygies-1",
    "page": "Finitely generated modules",
    "title": "Syzygies",
    "category": "section",
    "text": "syz(::smodule)ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nv1 = vector(R, (x + 1)*y, (x*y + 1)*y, y)\nv2 = vector(R, (x + 1)*x, (x*y + 1)*x, x)\n\nM = Singular.Module(R, v1, v2)\n\nZ = syz(M)"
},

{
    "location": "module.html#Singular.sres-Union{Tuple{Singular.smodule{T},Int64}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Finitely generated modules",
    "title": "Singular.sres",
    "category": "Method",
    "text": "sres{T <: Nemo.RingElem}(I::smodule{T}, max_length::Int)\n\nCompute a free resolution of the given module I of length up to the given maximum length. If max_length is set to zero, a full length free resolution is computed. Each element of the resolution is itself a module.\n\n\n\n"
},

{
    "location": "module.html#Free-resolutions-1",
    "page": "Finitely generated modules",
    "title": "Free resolutions",
    "category": "section",
    "text": "sres{T <: Nemo.RingElem}(::smodule{T}, ::Int)ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nv1 = vector(R, x + 1, x*y + 1, y)\nv2 = vector(R, x^2 + 1, 2x + 3y, x)\n\nM = std(Singular.Module(R, v1, v2))\n\nF = sres(M, 0)\n\nM1 = Singular.Matrix(M)\nM2 = Singular.Matrix(F[2])\n\n# test we have a complex\niszero(M1*M2)"
},

{
    "location": "vector.html#",
    "page": "Free modules and vectors",
    "title": "Free modules and vectors",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "vector.html#Free-modules-and-vectors-1",
    "page": "Free modules and vectors",
    "title": "Free modules and vectors",
    "category": "section",
    "text": "As generators of finitely generated modules in Singular.jl are given as submodule of free modules over a polynomial ring R, Singular.jl supports creation of the free module R^n and vectors of length n in such a module.The Singular.jl type for a vector is svector{T}. For the most part, these exist to help interact with the smodule{T} type provided by Singular.The types of vectors and associated parent objects are given in the following table according to the library provding them.Library Element type Parent type\nSingular svector{T} Singular.FreeMod{T}These types are parameterised by the type of elements in the polynomial ring R.All free module types belong directly to the abstract type Module{T} and vector types belong directly to ModuleElem{T}."
},

{
    "location": "vector.html#Free-module-and-vector-functionality-1",
    "page": "Free modules and vectors",
    "title": "Free module and vector functionality",
    "category": "section",
    "text": "Singular.jl modules implement standard operations one would expect on vectors and their associated parent modules.These include:Operations common to all AbstractAlgebra objects, such as parent, base_ring, elem_type, parent_type, parent, deepcopy, etc.\nArithmetic operations on vectors: (unary) -, +, -\nScalar multiplication of vectors by polynomials, constants and integers\nComparison operators: ==Below, we describe all of the functionality for Singular.jl free modules that is not included in this list of basic operations."
},

{
    "location": "vector.html#Constructors-1",
    "page": "Free modules and vectors",
    "title": "Constructors",
    "category": "section",
    "text": "Given a Singular polynomial ring R and a rank n, the following constructors are available for creating free modules.FreeModule{T <: AbstractAlgebra.RingElem}(R::PolyRing{T}, n::Int)Construct the free module R^n over the polynomial ring R. Elements of the module returned by this function are vectors of length n, with entries in R.vector{T <: AbstractAlgebra.RingElem}(R::PolyRing{T}, coords::spoly{T}...)Construct the vector whose coordinates (which are elements of R) are the given parameters.ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nM = FreeModule(R, 3)\nv2 = M([x + 1, x*y + 1, y])\n\nv1 = vector(R, x + 1, x*y + 1, y)"
},

{
    "location": "vector.html#Base.LinAlg.rank-Tuple{Singular.FreeMod}",
    "page": "Free modules and vectors",
    "title": "Base.LinAlg.rank",
    "category": "Method",
    "text": "rank(M::FreeMod)\n\nReturn the rank of the given free module.\n\n\n\n"
},

{
    "location": "vector.html#AbstractAlgebra.Generic.gens-Union{Tuple{Singular.FreeMod{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Free modules and vectors",
    "title": "AbstractAlgebra.Generic.gens",
    "category": "Method",
    "text": "gens{T <: AbstractAlgebra.RingElem}(M::FreeMod{T})\n\nReturn a Julia array whose entries are the generators of the given free module.\n\n\n\n"
},

{
    "location": "vector.html#Basic-manipulation-1",
    "page": "Free modules and vectors",
    "title": "Basic manipulation",
    "category": "section",
    "text": "rank(::FreeMod)gens{T <: AbstractAlgebra.RingElem}(::FreeMod{T})ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nM = FreeModule(R, 5)\n\nv = gens(M)\nr = rank(M)"
},

{
    "location": "vector.html#Conversions-1",
    "page": "Free modules and vectors",
    "title": "Conversions",
    "category": "section",
    "text": "To convert the internal Singular representation of an svector{T} to a Julia array whose entries have type T, we have the following conversion routine.Array{T <: Nemo.RingElem}(v::svector{spoly{T}})ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nv1 = vector(R, x + 1, x*y + 1, y)\n\nV = Array(v1)"
},

{
    "location": "resolution.html#",
    "page": "Resolutions",
    "title": "Resolutions",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "resolution.html#Resolutions-1",
    "page": "Resolutions",
    "title": "Resolutions",
    "category": "section",
    "text": "Functions for creating free resolutions of modules and ideals in Singular.jl return a special Singular object of type sresolution{T}. The support in Singular.jl for this type primarily exists to allow interaction with such resolutions. Free resolutions can have the property of being minimal, which is specified by the minimal field of the sresolution{T} type.Resolution objects have a parent object which represents the set of resolutions they belong to, the data for which is given by the polynomial ring R over which the modules in the resolution are defined.The types of resolutions and associated parent objects are given in the following table according to the library provding them.Library Element type Parent type\nSingular sresolution{T} Singular.ResolutionSet{T}These types are parameterised by the type of elements in the polynomial ring R over which the modules belonging to the resolution are defined.All resolution types belong directly to the abstract type SetElem and all the resolution set parent object types belong to the abstract type Set."
},

{
    "location": "resolution.html#Resolution-functionality-1",
    "page": "Resolutions",
    "title": "Resolution functionality",
    "category": "section",
    "text": "Singular.jl resolutions implement standard operations one would expect on all AbstractAlgebra compatible objects. These include:Operations common to all AbstractAlgebra objects, such as parent, base_ring, elem_type, parent_type, parent, deepcopy, etc.Below, we describe all of the functionality for Singular.jl resolutions that is not included in this list of basic operations."
},

{
    "location": "resolution.html#Constructors-1",
    "page": "Resolutions",
    "title": "Constructors",
    "category": "section",
    "text": "Resolutions can currently only be created by taking the free resolution of an ideal or module over a polynomial ring, as described in the relevant sections of the documentation.Alternatively, resolutions can be refined to minimal resolutions, as described below.Other than this, there are currently no additional ways to create resolutions in Singular.jl."
},

{
    "location": "resolution.html#Base.length-Tuple{Singular.sresolution}",
    "page": "Resolutions",
    "title": "Base.length",
    "category": "Method",
    "text": "length(r::sresolution)\n\nReturn the length of the resolution. This is what is mathematically meant by the length of a resolution. Over a field, this should be at most the number of variables in the polynomial ring.\n\n\n\n"
},

{
    "location": "resolution.html#Basic-manipulation-1",
    "page": "Resolutions",
    "title": "Basic manipulation",
    "category": "section",
    "text": "length(::sresolution)Singular.jl overloads the getindex function so that one can access the modules in a resolution F.F[n::Int]ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x*y + 1, x^2 + 1)\nF = fres(std(I), 0)\n\nn = length(F)\nM1 = F[1]"
},

{
    "location": "resolution.html#Singular.betti-Tuple{Singular.sresolution}",
    "page": "Resolutions",
    "title": "Singular.betti",
    "category": "Method",
    "text": "betti(r::sresolution)\n\nReturn the Betti numbers, i.e. the ranks of the free modules in the given free resolution. These are returned as a Julia array of Ints.\n\n\n\n"
},

{
    "location": "resolution.html#Betti-numbers-1",
    "page": "Resolutions",
    "title": "Betti numbers",
    "category": "section",
    "text": "betti(::sresolution)ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x*y + 1, x^2 + 1)\nF = fres(std(I), 3)\nM = minres(F)\n\nB = betti(M)"
},

{
    "location": "resolution.html#Singular.minres-Union{Tuple{Singular.sresolution{T}}, Tuple{T}} where T<:AbstractAlgebra.RingElem",
    "page": "Resolutions",
    "title": "Singular.minres",
    "category": "Method",
    "text": "minres{T <: AbstractAlgebra.RingElem}(r::sresolution{T})\n\nReturn a minimal free resolution, given any free resolution. If the supplied resolution is already minimal, it may be returned without making a copy.\n\n\n\n"
},

{
    "location": "resolution.html#Minimal-resolutions-1",
    "page": "Resolutions",
    "title": "Minimal resolutions",
    "category": "section",
    "text": "minres{T <: AbstractAlgebra.RingElem}(::sresolution{T})ExamplesR, (x, y) = PolynomialRing(QQ, [\"x\", \"y\"])\n\nI = Ideal(R, x*y + 1, x^2 + 1)\nF = fres(std(I), 3)\nM = minres(F)"
},

{
    "location": "matrix.html#",
    "page": "Matrices",
    "title": "Matrices",
    "category": "page",
    "text": "CurrentModule = Singular"
},

{
    "location": "matrix.html#Matrices-1",
    "page": "Matrices",
    "title": "Matrices",
    "category": "section",
    "text": "Singular internally allows for matrices over polynomial rings to be created extremely efficiently from ideals and modules (often without copying data). This allows for introspection of modules and operations that can be expressed in terms of matrices (e.g. composition of R-module homomorphisms) to be computed, at a low level.The default matrix type in Singular.jl is the smatrix type.Matrix objects have a parent object which represents the space of matrices they belong to, the data for which is given by the polynomial ring R over which the matrices are defined, and the number of rows and columns of the matrices in the space.The types of matrices and associated parent objects are given in the following table according to the library provding them.Library Element type Parent type\nSingular smatrix{T} Singular.MatrixSpace{T}These types are parameterised by the type of elements in the polynomial ring R over which the matrices are defined.All matrix types belong directly to the abstract type SetElem and all the matrix space parent object types belong to the abstract type Set."
},

{
    "location": "matrix.html#Matrix-functionality-1",
    "page": "Matrices",
    "title": "Matrix functionality",
    "category": "section",
    "text": "Singular.jl matrices implement standard operations one would expect. These include:Operations common to all AbstractAlgebra objects, such as parent, base_ring, elem_type, parent_type, parent, deepcopy, etc.The following parts of the Matrix interface from AbstractAlgebra are also implemented:matrix dimensions: nrows, ncols\naccessing elements of a matrix M[i, j]\niszero\narithmetic operations: +, -, *\ncomparison: ==No other functionality is provided at this time."
},

]}
