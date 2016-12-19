###############################################################################
#
#   Coefficients
#
###############################################################################

typealias coeffs pcpp"n_Procs_s"

typealias n_coeffType Cxx.CppEnum{:n_coeffType}

typealias number pcpp"snumber"

typealias number_ref Ref{number}

typealias __mpz_struct pcpp"__mpz_struct"

typealias mpz_t pcpp"mpz_t"

###############################################################################
#
#   Polynomials/vectors
#
###############################################################################

typealias ring pcpp"ip_sring"

typealias poly pcpp"spolyrec"

typealias poly_ref Ref{poly}

typealias vector pcpp"spolyrec"

typealias rRingOrder_t Cxx.CppEnum{:rRingOrder_t, Cuint}

###############################################################################
#
#   Ideals
#
###############################################################################

typealias ideal pcpp"sip_sideal"

typealias ideal_ref Ref{ideal}

###############################################################################
#
#   Matrices
#
###############################################################################

typealias matrix pcpp"ip_smatrix"

###############################################################################
#
#   Resolvente (module list)
#
###############################################################################

typealias resolvente Cxx.CppPtr{pcpp"sip_sideal", (false, false, false)}

###############################################################################
#
#   Syzygies
#
###############################################################################

typealias syStrategy pcpp"ssyStrategy" 
