#include "includes.h"
#include "coeffs.h"
#include "rings.h"

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

  singular_define_coeffs(Singular);
  singular_define_rings(Singular);

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
