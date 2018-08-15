#include <string>
#include <cstring>
#include <iostream>
#include <tuple>

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

auto id_fres_helper(sip_sideal* I, int n, std::string method, ring R){
    auto origin = currRing;
    rChangeCurrRing(R);
    syStrategy s = syFrank(I, n, method.c_str());
    rChangeCurrRing(origin);
    auto r = s->minres;
    bool minimal = true;
    if(r==NULL){
        r = s->fullres;
        minimal = false;
    }
    return std::make_tuple(reinterpret_cast<void*>(r),s->length,minimal);
}

auto rDefault_helper(coeffs cf, jlcxx::ArrayRef<std::string> vars, rRingOrder_t ord){
    auto len = vars.size();
    char** vars_ptr  = new char*[len];
    for(int i = 0;i<len; i++){
        vars_ptr[i] = new char[vars[i].length()+1];
        std::strcpy(vars_ptr[i],vars[i].c_str());
    }
    auto r = rDefault(cf,len,vars_ptr,ord);
    delete[] vars_ptr;
    r->ShortOut = 0;
    return r;
}

// typedef snumber* snumberptr;

JULIA_CPP_MODULE_BEGIN(registry)
  jlcxx::Module& Singular = registry.create_module("libSingular");

  Singular.add_type<n_Procs_s>("coeffs");
  Singular.add_bits<n_coeffType>("n_coeffType");
  Singular.set_const("n_Z",n_Z);
  Singular.set_const("n_Q",n_Q);
  Singular.set_const("n_Zn",n_Zn);
  Singular.set_const("n_Zp",n_Zp);
  Singular.set_const("n_GF",n_GF);
  Singular.set_const("n_unknown",n_unknown);
  Singular.add_type<snumber>("number");
  Singular.add_type<__mpz_struct>("__mpz_struct");
  Singular.add_type<ip_sring>("ring");
  Singular.add_type<spolyrec>("poly");
  // Singular.add_type<nMapFunc>("nMapFunc");
  // Singular.add_type<spolyrec>("vector");
  Singular.add_bits<rRingOrder_t>("rRingOrder_t");
  Singular.add_type<sip_sideal>("ideal");
  Singular.add_type<ip_smatrix>("ip_smatrix");
  Singular.add_type<ssyStrategy>("syStrategy");

  /* monomial orderings */
  Singular.set_const("ringorder_no", ringorder_no);
  Singular.set_const("ringorder_lp", ringorder_lp);
  Singular.set_const("ringorder_rp", ringorder_rp);
  Singular.set_const("ringorder_dp", ringorder_dp);
  Singular.set_const("ringorder_Dp", ringorder_Dp);
  Singular.set_const("ringorder_ls", ringorder_ls);
  Singular.set_const("ringorder_rs", ringorder_rs);
  Singular.set_const("ringorder_ds", ringorder_ds);
  Singular.set_const("ringorder_Ds", ringorder_Ds);
  Singular.set_const("ringorder_c", ringorder_c);
  Singular.set_const("ringorder_C", ringorder_C);

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
  Singular.method("n_SetMap",[]( const coeffs x, const coeffs y){ return reinterpret_cast<void*>(n_SetMap(x,y)); });

  Singular.method("nApplyMapFunc",[]( void* map, snumber* x, coeffs a, coeffs b ){ return reinterpret_cast<nMapFunc>(map)(x,a,b); });

  Singular.method("n_Init",[]( long x, coeffs n){ return n_Init(x,n); });

  Singular.method("n_Copy",[]( snumber* x, const coeffs n){ return n_Copy(x,n); });

  Singular.method("nCoeff_has_simple_Alloc",[]( coeffs x ){ return nCoeff_has_simple_Alloc( x ) > 0; });

  Singular.method("n_InitMPZ_internal",&n_InitMPZ);

  Singular.method("n_Delete",[]( snumber* n, coeffs cf) { number* t = &n; if( n != NULL ){ n_Delete(t, cf); } });

  Singular.method("n_Delete_Q",[]( void* n, coeffs cf ) { number tt = reinterpret_cast<number>(n); number* t = &tt; if( n != NULL ){ n_Delete(t, cf); } });

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

  Singular.method("StringEndS",[](){ char* m = StringEndS(); std::string s(m); omFree(m); return s; });

  Singular.method("omAlloc0",[]( size_t siz ){ return (void*) omAlloc0(siz); });

  Singular.method("omFree_internal",[]( void* m ){ omFree(m); });

  /* Setting a Ptr{number} to a number */

  Singular.method("setindex_internal",[]( void* x, snumber* y ){ *reinterpret_cast<snumber**>(x) = y; });

  Singular.method("id_fres",&id_fres_helper);
  
  /****************************
   ** from coeffs.jl
   ***************************/
  Singular.method("rDefault",&rDefault_helper);
  Singular.method("rDelete",&rDelete);
  Singular.method("rString",[](ip_sring* r){auto s = rString(r); return std::string(s);});
  Singular.method("rChar",&rChar);
  Singular.method("rGetVar",&rGetVar);
  Singular.method("rVar",&rVar);
  Singular.method("rGetExpSize",[]( unsigned long bitmask, int N){ int bits; return static_cast<unsigned int>(rGetExpSize(bitmask,bits,N));});
  Singular.method("rHasGlobalOrdering",&rHasGlobalOrdering);

  /****************************
   ** from ideals.jl
   ***************************/

  Singular.method("id_Delete",[]( ideal m, ring n ) { return id_Delete(&m,n); });

  Singular.method("id_Copy",&id_Copy);

  Singular.method("idInit",&idInit);

  Singular.method("setindex!",[]( ideal r, poly n, int o) { r->m[o] = n; });

  Singular.method("getindex",[]( ideal r, int o) { (poly) (r->m[o]); });
  
  Singular.method("idIs0",&idIs0);

  Singular.method("id_IsConstant",&id_IsConstant);

  Singular.method("id_IsZeroDim",&id_IsZeroDim);

  Singular.method("idElem",&idElem);

  Singular.method("id_Normalize",&id_Normalize);

  Singular.method("id_Head",&id_Head);

  Singular.method("id_MaxIdeal",[](int m, ring n) { return id_MaxIdeal(m,n); });

  Singular.method("id_Add",&id_Add);

  Singular.method("id_Mult",&id_Mult);

  Singular.method("id_Power",&id_Power);

  Singular.method("id_IsEqual",[](ideal m, ideal n, ring o) { mp_Equal((ip_smatrix *) m, (ip_smatrix *) n, o); });

  Singular.method("id_FreeModule",&id_FreeModule);

  Singular.method("idSkipZeroes",&idSkipZeroes);

  Singular.method("ngens",[](ideal m) { (int) IDELEMS(m); });

  Singular.method("rank",[](ideal m ) { (int) m->rank; });

  Singular.method("id_Quotient",[](ideal a, ideal b, bool c, ring d) { 
    const ring origin = currRing;
          rChangeCurrRing(d);
          ideal id = idQuot(a, b, c, TRUE);
          rChangeCurrRing(origin);
          id; });

  Singular.method("id_Intersection",[](ideal a, ideal b, ring c) { 
    const ring origin = currRing;
          rChangeCurrRing(c);
          ideal id = idSect(a, b);
          rChangeCurrRing(origin);
          id; });

 /* Singular.method("id_Slimgb",[](ideal a, ring b, bool c= false) {
    if (c)
      crbit = Sy_bit(OPT_REDSB);
   else
      crbit = unsigned int(0);
   ideal id = NULL;
          if (!idIs0(a))
          {
             intvec * n = NULL;
             tHomog h = testHomog;
             const ring origin = currRing;
             unsigned int save_opt = si_opt_1;
             si_opt_1 |= crbit;
             rChangeCurrRing(b);
             id = t_rep_gb(b, a, a->rank);
             si_opt_1 = save_opt;
             rChangeCurrRing(origin);
             if (n != NULL)
                delete n;
          } else
             id = idInit(0, a->rank);
          id;
       }); 
*/

  Singular.method("id_Syzygies",[](ideal m, ring o) { 
          ideal id = NULL;
          intvec * n = NULL;
          tHomog h = testHomog;
          const ring origin = currRing;
          rChangeCurrRing(o);
          id = idSyzygies(m, h, &n);
          rChangeCurrRing(origin);
          if (n != NULL)
             delete n;
          id; });
/*
  Singular.method("id_sres",[](ideal m, int n, String s, ring o) {
                s = const ring origin = currRing;
         rChangeCurrRing(o);
         syStrategy s = sySchreyer(m, n);
         rChangeCurrRing(origin);
         s;
   r = s->minres;
   minimal = true
   if r == C_NULL
      r = s->fullres;
      minimal = false
   end
   length = s->length;
   return r, Int(length), minimal;

});
*/

  Singular.method("id_Eliminate",[](ideal m, poly p, ring o) {
          const ring origin = currRing;
          rChangeCurrRing(o);
          ideal res = idElimination(m, p);
          rChangeCurrRing(origin);
          res; });

JULIA_CPP_MODULE_END
