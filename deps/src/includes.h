#ifndef MAIN_INCLUDE_FILE
#define MAIN_INCLUDE_FILE

#include <string>
#include <cstring>
#include <iostream>
#include <tuple>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/const_array.hpp"
#include "jlcxx/array.hpp"
#include "jlcxx/tuple.hpp"

#include <gmp.h>
#include <omalloc/omalloc.h>
#include <misc/intvec.h>
#include <misc/auxiliary.h>
#include <reporter/reporter.h>
// #include <feFopen.h>
#include <coeffs/coeffs.h>
#include <coeffs/longrat.h>
#include <polys/clapsing.h>
#include <coeffs/bigintmat.h>
#include <coeffs/rmodulon.h>
#include <polys/monomials/ring.h>
#include <polys/monomials/p_polys.h>
#include <polys/simpleideals.h>
#include <kernel/GBEngine/kstd1.h>
#include <kernel/GBEngine/syz.h>
#include <kernel/GBEngine/tgb.h>
#include <kernel/ideals.h>
#include <kernel/polys.h>
#include <kernel/combinatorics/stairc.h>
#include <Singular/grammar.h>
#include <Singular/libsingular.h>
#include <Singular/fevoices.h>
#include <Singular/ipshell.h>
#include <Singular/ipid.h>
#include <Singular/subexpr.h>
#include <Singular/lists.h>
#include <Singular/idrec.h>
#include <Singular/tok.h>
#include <Singular/links/silink.h>
#include <Singular/fehelp.h>

namespace jlcxx {
template <> struct IsBits<n_coeffType> : std::true_type {
};
template <> struct IsBits<rRingOrder_t> : std::true_type {
};
}    // namespace jlcxx

#endif
