#include <string>
#include <iostream>

#include "jlcxx/jlcxx.hpp"

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
    template<> struct IsBits<n_coeffType> : std::true_type {};
    template<> struct IsBits<rRingOrder_t> : std::true_type {};
}

JULIA_CPP_MODULE_BEGIN(registry)
  jlcxx::Module& Singular = registry.create_module("libSingular");

  Singular.add_type<n_Procs_s>("coeffs");
  Singular.add_bits<n_coeffType>("n_coeffType");
  Singular.set_const("n_Z",n_Z);
  Singular.add_type<snumber>("number");
  Singular.add_type<__mpz_struct>("__mpz_struct");
  Singular.add_type<ip_sring>("ring");
  Singular.add_type<spolyrec>("poly");
  // Singular.add_type<spolyrec>("vector");
  Singular.add_bits<rRingOrder_t>("r_RingOrder_t");
  Singular.add_type<sip_sideal>("ideal");
  Singular.add_type<ip_smatrix>("ip_smatrix");
  Singular.add_type<ssyStrategy>("syStrategy");

  Singular.method("siInit",[](const char* path){ siInit(const_cast<char*>(path)); });
  Singular.method("versionString",[](){ return const_cast<const char*>(versionString()); });

  // from coeffs.jl
  Singular.method("nInitChar",&nInitChar);
  Singular.method("n_GetChar",[]( const coeffs n ){ return n_GetChar(n);});
  Singular.method("nCopyCoeff",&nCopyCoeff);
  Singular.method("nKillChar",&nKillChar);
  // Singular.method("nApplyMapFunc",[]())
  // Singular.method("n_Init",[](long x, coeffs n){ n_Init);



JULIA_CPP_MODULE_END
