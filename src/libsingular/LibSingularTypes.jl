###############################################################################
#
#   Coefficients
#
###############################################################################

const number_ref = CxxRef{number}

const mpz_t = Ptr{__mpz_struct}

###############################################################################
#
#   Polynomials/vectors
#
###############################################################################

const poly_ref = CxxRef{poly}

const vector = poly

###############################################################################
#
#   Ideals
#
###############################################################################

const ideal_ref = CxxRef{ideal}

###############################################################################
#
#   Matrices
#
###############################################################################

const matrix_ref = CxxRef{ip_smatrix}

###############################################################################
#
#   Resolvente (module list)
#
###############################################################################

const resolvente = Ptr{ideal}

const syStrategy_ref = CxxRef{syStrategy}

###############################################################################
#
#   Syzygies
#
###############################################################################

###############################################################################
#
#   CoeffRingData
#
###############################################################################

mutable struct singular_coeff_ring_struct
    has_simple_alloc::Int64
    has_simple_inverse::Int64
    is_field::Int64
    is_domain::Int64
    ch::Int64
    data::Ptr{Cvoid}
    cfInit::Ptr{Cvoid}
    cfInt::Ptr{Cvoid}
    cfMPZ::Ptr{Cvoid}
    cfInpNeg::Ptr{Cvoid}
    cfCopy::Ptr{Cvoid}
    cfDelete::Ptr{Cvoid}
    cfAdd::Ptr{Cvoid}
    cfInpAdd::Ptr{Cvoid}
    cfSub::Ptr{Cvoid}
    cfMult::Ptr{Cvoid}
    cfInpMult::Ptr{Cvoid}
    cfDiv::Ptr{Cvoid}
    cfDivBy::Ptr{Cvoid}
    cfInvers::Ptr{Cvoid}
    cfGcd::Ptr{Cvoid}
    cfSubringGcd::Ptr{Cvoid}
    cfExtGcd::Ptr{Cvoid}
    cfGreater::Ptr{Cvoid}
    cfEqual::Ptr{Cvoid}
    cfIsZero::Ptr{Cvoid}
    cfIsOne::Ptr{Cvoid}
    cfIsMOne::Ptr{Cvoid}
    cfGreaterZero::Ptr{Cvoid}
    cfWriteLong::Ptr{Cvoid}
    cfCoeffWrite::Ptr{Cvoid}
end

function singular_coeff_ring_struct()
    singular_coeff_ring_struct(Tuple(cat([0 for i in 1:5], [Ptr{Cvoid}(0) for i in 1:26], dims = 1))...)
end
