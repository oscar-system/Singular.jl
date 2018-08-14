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

// typedef snumber* snumberptr;

JULIA_CPP_MODULE_BEGIN(registry)
  jlcxx::Module& Singular = registry.create_module("libSingular");

  Singular.add_type<n_Procs_s>("coeffs");
  Singular.add_bits<n_coeffType>("n_coeffType");
  Singular.set_const("n_Z",n_Z);
  Singular.set_const("n_Q",n_Q);
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

  /****************************
   ** from coeffs.jl
   ***************************/

  /* initialise a coefficient ring */
  Singular.method("nInitChar",&nInitChar);

  /* get the characteristic of a coefficient domain */
  Singular.method("n_GetChar",[]( coeffs n ){ return n_GetChar(n);});

  /* make a copy of a coefficient domain (actually just increments a reference count) */
  Singular.method("nCopyCoeff",&nCopyCoeff);

  /* kill a coefficient ring */
  Singular.method("nKillChar",&nKillChar);

  /* return a function to convert between rings */
  // Singular.method("n_SetMap",[]( const coeffs x, const coeffs y){ return n_SetMap(x,y); });

  // Singular.method("nApplyMapFunc",[]())

  Singular.method("n_Init",[]( long x, coeffs n){ return n_Init(x,n); });

  Singular.method("n_Copy",[]( snumber* x, const coeffs n){ return n_Copy(x,n); });

  Singular.method("nCoeff_has_simple_Alloc",[]( coeffs x ){ return nCoeff_has_simple_Alloc( x ) > 0; });

  Singular.method("n_InitMPZ_internal",&n_InitMPZ);

  Singular.method("n_Delete",[]( snumber n, coeffs cf) { number t = &n; if (t != NULL) n_Delete(&t, cf); });

  Singular.method("n_Write_internal",[]( snumber* x, coeffs cf, const int d ){ return n_Write(x,cf,d); });

  Singular.method("n_Add",[]( snumber* a, snumber* b, coeffs c ){ return n_Add(a,b,c); });

  Singular.method("n_Sub",[]( snumber* a, snumber* b, coeffs c ){ return n_Sub(a,b,c); });

  Singular.method("n_Mult",[]( snumber* a, snumber* b, coeffs c ){ return n_Mult(a,b,c); });

  Singular.method("n_Neg",[]( snumber* a, coeffs c ){ number nn = n_Copy(a, c); nn = n_InpNeg(nn, c); return nn; });

  Singular.method("n_Invers",[]( snumber* a, coeffs c ){ return n_Invers(a,c); });

  Singular.method("n_ExactDiv",[]( snumber* a, snumber* b, coeffs c ){ return n_ExactDiv(a,b,c); });

  Singular.method("n_Div",[]( snumber* a, snumber* b, coeffs c ){ number z = n_Div(a, b, c); n_Normalize(z, c); return z; });

  Singular.method("n_GetNumerator",[]( snumber* a, coeffs c){ return n_GetNumerator(a,c); });

  Singular.method("n_GetDenom",[]( snumber* a, coeffs c){ return n_GetDenom(a,c); });

  Singular.method("n_Power",[]( snumber* a, int b, coeffs c ){ number res; n_Power(a, b, &res, c); return res; });

  Singular.method("n_Gcd",[]( snumber* a, snumber* b, coeffs c ){ return n_Gcd(a,b,c); });

  Singular.method("n_SubringGcd",[]( snumber* a, snumber* b, coeffs c ){ return n_SubringGcd(a,b,c); });

  Singular.method("n_Lcm",[]( snumber* a, snumber* b, coeffs c ){ return n_Lcm(a,b,c); });

  Singular.method("n_ExtGcd_internal",[]( snumber* a, snumber* b, void* s, void* t, coeffs c ){
                                 return n_ExtGcd(a,b,
                                          reinterpret_cast<snumber**>(s),
                                          reinterpret_cast<snumber**>(t),c); });

  Singular.method("n_IsZero",[]( snumber* x, const coeffs n ){ return n_IsZero( x,n ) > 0; });

  Singular.method("n_IsOne",[]( snumber* x, const coeffs n ){ return n_IsOne( x,n ) > 0; });

  Singular.method("n_Greater",[]( snumber* x, snumber* y, const coeffs n ){ return n_Greater( x,y,n ) > 0; });

  Singular.method("n_Equal",[]( snumber* x, snumber* y, const coeffs n ){ return n_Equal( x,y,n ) > 0; });

  Singular.method("n_InpAdd",[]( snumber* x, snumber* y, const coeffs n ){ return n_InpAdd( x,y,n ); });

  Singular.method("n_InpMult",[]( snumber* x, snumber* y, const coeffs n ){ return n_InpMult( x,y,n ); });

  Singular.method("n_QuotRem_internal",[]( snumber* x, snumber* y, void* p, const coeffs n ){ return n_QuotRem( x,y,reinterpret_cast<snumber**>(p),n ); });

  Singular.method("n_Rem",[]( snumber* x, snumber* y, const coeffs n ){ number qq; return n_QuotRem( x,y, &qq, n ); });

  Singular.method("n_IntMod",&n_IntMod);

  Singular.method("n_Farey",&n_Farey);

  Singular.method("n_ChineseRemainderSym_internal",[]( void* x, void* y, int n, int sign_flag, coeffs c ){ CFArray inv_cache(n); return n_ChineseRemainderSym( reinterpret_cast<snumber**>(x),reinterpret_cast<snumber**>(y),n,sign_flag,inv_cache,c ); });

  Singular.method("n_Param",[]( int x, const coeffs n ){ return n_Param(x,n); });

  Singular.method("StringSetS_internal",[]( std::string m ){ return StringSetS(m.c_str()); });

  Singular.method("StringEndS",[](){ return std::string(StringEndS()); });

  Singular.method("omAlloc0",[]( size_t siz ){ return (void*) omAlloc0(siz); });

  Singular.method("omFree_internal",[]( void* m ){ omFree(m); });

  /* Setting a Ptr{number} to a number */

  Singular.method("setindex_internal",[]( void* x, snumber* y ){ *reinterpret_cast<snumber**>(x) = y; });

JULIA_CPP_MODULE_END
