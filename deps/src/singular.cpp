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



auto id_sres_helper(sip_sideal* m, int n, ring R) {
  auto origin = currRing;
  rChangeCurrRing(R);
  syStrategy s = sySchreyer(m, n);
  rChangeCurrRing(origin);
  auto r = s->minres;
  bool minimal = true;
  if (r == NULL){
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

ideal id_Syzygies_internal(ideal m, ring o) 
{ 
  ideal id = NULL;
  intvec * n = NULL;
  tHomog h = testHomog;
  const ring origin = currRing;
  rChangeCurrRing(o);
  id = idSyzygies(m, h, &n);
  rChangeCurrRing(origin);
  if (n != NULL)
     delete n;
  return id; 
}

auto id_Slimgb_helper(ideal a, ring b,bool complete_reduction) {
//  bool complete_reduction= false;
  unsigned int crbit;
  if (complete_reduction == false)
	auto crbit = Sy_bit(OPT_REDSB);
  else
	crbit = 0;
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
          return id;
} 

auto id_Std_helper(ideal a, ring b) {
  bool complete_reduction= false;
  unsigned int crbit;
  if (complete_reduction == false)
      crbit = Sy_bit(OPT_REDSB);
  else
	crbit = 0;
   	ideal id = NULL;
          if (!idIs0(a))
          {
             intvec * n = NULL;
             tHomog h = testHomog;
             const ring origin = currRing;
             unsigned int save_opt = si_opt_1;
             si_opt_1 |= crbit;
             rChangeCurrRing(b);
             id = kStd(a, b->qideal, h, &n);
             si_opt_1 = save_opt;
             rChangeCurrRing(origin);
             if (n != NULL)
                delete n;
          } else
             id = idInit(0, a->rank);
          return id;
}



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
  Singular.method("rBitmask",[](ip_sring* r){ return (unsigned int)r->bitmask; });
  Singular.method("p_Delete",[](spolyrec* p, ip_sring* r){ return p_Delete(&p,r);});
  Singular.method("p_Copy",[](spolyrec* p, ip_sring* r){ return p_Copy(p,r);});
  Singular.method("p_IsOne",[](spolyrec* p, ip_sring* r){ return p_IsOne(p,r);});
  Singular.method("p_One",[](ip_sring* r){ return p_One(r);});
  Singular.method("p_Init",[](ip_sring*r ){ return p_Init(r);} );
  Singular.method("p_IsUnit",[](spolyrec* p, ip_sring* r){ return p_IsUnit(p,r);});
  Singular.method("p_GetExp",[](spolyrec* p, int i, ip_sring* r){ return p_GetExp(p,i,r);});
  Singular.method("p_GetComp",[](spolyrec* p, ip_sring* r){ return p_GetComp(p,r);});
  Singular.method("p_String",[](spolyrec* p, ip_sring* r){ std::string s(p_String(p,r)); return s; });
  Singular.method("p_ISet",[](long i, ip_sring* r){ return p_ISet(i,r);});
  Singular.method("p_NSet",[](snumber* p, ip_sring* r){ return p_NSet(p,r);});
  Singular.method("pLength",[](spolyrec* p){ return pLength(p);});
  Singular.method("pNext",[](spolyrec* p){ return pNext(p);});
  Singular.method("p_Neg",[](spolyrec* p, ip_sring* r){ return p_Neg(p,r); });
  Singular.method("pGetCoeff",[](spolyrec* p){ return pGetCoeff(p);});
  Singular.method("pSetCoeff",[](spolyrec* p, long c, ip_sring* r){ number n = n_Init(c,r); return p_SetCoeff(p,n,r); });
  Singular.method("pSetCoeff0",[](spolyrec* p, long c, ip_sring* r){ number n = n_Init(c,r); return p_SetCoeff0(p,n,r); });
  Singular.method("pLDeg",[](spolyrec* a, ip_sring* r){ long res; int dummy; if(a==NULL){ res = r->pLDeg(a,&dummy,r);}else{ res = -1;} return res; });
  Singular.method("p_Add_q",[](spolyrec* p, spolyrec* q, ip_sring* r){ return p_Add_q(p,q,r); });
  Singular.method("p_Sub",[](spolyrec* p, spolyrec* q, ip_sring* r){ return p_Sub(p,q,r); });
  Singular.method("p_Mult_q",[](spolyrec* p, spolyrec* q, ip_sring* r){ return p_Mult_q(p,q,r); });
  Singular.method("p_Power",[](spolyrec* p, int q, ip_sring* r){ return p_Power(p,q,r); });
  Singular.method("p_EqualPolys",[](spolyrec* p, spolyrec* q, ip_sring* r){ return p_EqualPolys(p,q,r); });
  Singular.method("p_Divide",[](spolyrec* p, spolyrec* q, ip_sring* r){ return p_Divide(p,q,r); });
  Singular.method("singclap_gcd",[](spolyrec* p, spolyrec* q, ip_sring* r){ return singclap_gcd(p,q,r); });
  Singular.method("singclap_extgcd",[]( spolyrec* a, spolyrec* b, spolyrec* res, spolyrec* s, spolyrec* t, ip_sring* r ){
                                 return singclap_extgcd(a,b,
                                          res,
                                          s,
                                          t,r); });
  Singular.method("p_Content",[](spolyrec* p, ip_sring* r){ return p_Content(p,r); });
  Singular.method("p_GetExpVL_internal",[](spolyrec* p, long* ev, ip_sring* r){return p_GetExpVL(p,ev,r);});
  Singular.method("p_SetExpV_internal",[](spolyrec* p, int* ev, ip_sring* r){return p_SetExpV(p,ev,r);});
  Singular.method("p_Reduce",[](spolyrec* p, sip_sideal* G, ip_sring* R){
          const ring origin = currRing;
          rChangeCurrRing(R);
          poly res = kNF(G, R->qideal, p);
          rChangeCurrRing(origin);
          return res;
  });
  Singular.method("p_Reduce",[](sip_sideal* p, sip_sideal* G, ip_sring* R){
          const ring origin = currRing;
          rChangeCurrRing(R);
          ideal res = kNF(G, R->qideal, p);
          rChangeCurrRing(origin);
          return res;
  });

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

  Singular.method("id_Syzygies",&id_Syzygies_internal);

  Singular.method("id_fres",&id_fres_helper);

  Singular.method("id_sres",&id_sres_helper);

  Singular.method("id_Slimgb_helper",&id_Slimgb_helper);

  Singular.method("id_Std",&id_Std_helper);

  Singular.method("id_Eliminate",[](ideal m, poly p, ring o) {
          const ring origin = currRing;
          rChangeCurrRing(o);
          ideal res = idElimination(m, p);
          rChangeCurrRing(origin);
          res; });

  Singular.method("id_Satstd",&id_Satstd);

  Singular.method("id_Array2Vector",[](void* p, int a, ring o) { return id_Array2Vector(reinterpret_cast<spolyrec**>(p), a, o); });

  Singular.method("p_Vector2Array",[](poly p, void* a, int b, ring o) { return p_Vec2Array(p, reinterpret_cast<spolyrec**>(a), b, o); });







  /****************************
   ** from matrices.jl
   ***************************/

  Singular.method("ncols",[](matrix I){ return (int) MATCOLS(I);});

  Singular.method("nrows",[](matrix I){ return (int) MATROWS(I); });

  Singular.method("id_Module2Matrix",&id_Module2Matrix);

  Singular.method("getindex",[](matrix M, int i, int j){ return (poly) MATELEM(M, i, j); });

  Singular.method("mp_Copy",[](matrix M, ring R){ return mp_Copy(M, R);});

  Singular.method("mp_Delete",[](matrix M, ring R){ return mp_Delete(&M, R);});

  Singular.method("mp_Add",&mp_Add);

  Singular.method("mp_Sub",&mp_Sub);

  Singular.method("mp_Mult",&mp_Mult);

  Singular.method("mp_Equal",&mp_Equal);

 //Singular.method("iiStringMatrix",[](matrix I, int d, ring o){ return iiStringMatrix(I, d, o); });



  /****************************
   ** from resolutions.jl
   ***************************/

/*Singular.method("res_Delete",[](resolvente ra, int len, ring o){ 
	for (int i = 0; i < len; i++) {
          id_Delete(ra + i, o);
          omFreeSize((ADDRESS) ra, (len + 1)*sizeof(ideal));}});

*/

/*  Singular.method("res_Copy",[](resolvente ra, int len, ring o){ 
  resolvente res = (resolvente) omAlloc0((len + 1)*sizeof(ideal));
          rChangeCurrRing(o);
          for (int i = len - 1; i >= 0; i--)
          {
             if (ra[i] != NULL)
                res[i] = id_Copy(ra[i], o);
          }
          return res; });
*/


//  Singular.method("getindex",[](resolvente ra, int k){ return (ideal) ra[i];});

/*  Singular.method("syMinimize",[](resolvente ra, int len, ring o){
	  const ring origin = currRing;
          syStrategy temp = (syStrategy) omAlloc0(sizeof(ssyStrategy));
          resolvente result;
          rChangeCurrRing(o);
          temp->fullres = (resolvente) omAlloc0((len + 1)*sizeof(ideal));
          for (int i = len - 1; i >= 0; i--)
          {
             if (ra[i] != NULL)
                temp->fullres[i] = idCopy(ra[i]);
          }
          temp->length = len;
          syMinimize(temp);
          result = temp->minres;
          temp->minres = NULL;
          // syMinimize increments this as it returns a value we ignore 
          temp->references--;
          syKillComputation(temp, o);
          rChangeCurrRing(origin);
          return result;});



*/

/*
  Singular.method("syBetti",[](resolvente rs,int len, ring o){
   const ring origin = currRing;
   rChangeCurrRing(o);
   int dummy;
   intvec *iv = syBetti(rs, len, &dummy, NULL, FALSE, NULL);
   rChangeCurrRing(origin);
   return iv;
   int nrows = iv->rows();
   int ncols = iv->cols();
   auto betti = (int *)malloc(ncols*nrows*sizeof(int));
   	for (int i = 0; i < ncols; i++) {
            for (int j = 0; j < nrows; j++) {
               betti[i*nrows+j] = IMATELEM(*iv, j+1, i+1);
            }
         }
    delete(iv);
    return &betti[0];
  
   unsafe_wrap(Array, betti, (nrows, ncols), true);});
*/
JULIA_CPP_MODULE_END
