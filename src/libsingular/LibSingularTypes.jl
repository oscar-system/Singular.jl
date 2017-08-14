###############################################################################
#
#   Coefficients
#
###############################################################################

const coeffs = @pcpp_str "n_Procs_s"

const n_coeffType = Cxx.CppEnum{:n_coeffType}

const number = @pcpp_str "snumber"

const number_ref = Ref{number}

const __mpz_struct = @pcpp_str "__mpz_struct"

const mpz_t = @pcpp_str "mpz_t"

###############################################################################
#
#   Polynomials/vectors
#
###############################################################################

const ring = @pcpp_str "ip_sring"

const poly = @pcpp_str "spolyrec"

const poly_ref = Ref{poly}

const vector = @pcpp_str "spolyrec"

const rRingOrder_t = Cxx.CppEnum{:rRingOrder_t, Cuint}

###############################################################################
#
#   Ideals
#
###############################################################################

const ideal = @pcpp_str "sip_sideal"

const ideal_ref = Ref{ideal}

###############################################################################
#
#   Matrices
#
###############################################################################

const matrix = @pcpp_str "ip_smatrix"

###############################################################################
#
#   Resolvente (module list)
#
###############################################################################

const resolvente = Cxx.CppPtr{@pcpp_str "sip_sideal", (false, false, false)}

###############################################################################
#
#   Syzygies
#
###############################################################################

const syStrategy = @pcpp_str "ssyStrategy" 
