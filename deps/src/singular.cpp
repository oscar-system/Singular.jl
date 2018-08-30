#include "includes.h"
#include "coeffs.h"
#include "rings.h"
#include "ideals.h"
#include "matrices.h"



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
  Singular.add_type<sip_smap>("sip_smap");

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
  singular_define_ideals(Singular);
  singular_define_matrices(Singular);

  /****************************
   ** from resolutions.jl
   ***************************/

  Singular.method("res_Delete_helper",[](void* ra_void, int len, ring o){ 
    auto ra = reinterpret_cast<resolvente>(ra_void);
    for (int i = 0; i < len; i++) {
          id_Delete(ra + i, o);
          omFreeSize((ADDRESS) ra, (len + 1)*sizeof(ideal));}});

  Singular.method("res_Copy",[](void* ra_void, int len, ring o){ 
    auto ra = reinterpret_cast<resolvente>(ra_void);  
    resolvente res = (resolvente) omAlloc0((len + 1)*sizeof(ideal));
          rChangeCurrRing(o);
          for (int i = len - 1; i >= 0; i--)
          {
             if (ra[i] != NULL)
                res[i] = id_Copy(ra[i], o);
          }
          return res; });



  Singular.method("getindex",[](void* ra_void, int k){ 
    auto ra = reinterpret_cast<resolvente>(ra_void);
    return (ideal) ra[k];});

  Singular.method("syMinimize",[](void* ra_void, int len, ring o){
    auto ra = reinterpret_cast<resolvente>(ra_void); 
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



  Singular.method("syBetti",[](void* ra_void, int len, ring o){
    auto ra = reinterpret_cast<resolvente>(ra_void); 
    const ring origin = currRing;
    rChangeCurrRing(o);
    int dummy;
    intvec *iv = syBetti(ra, len, &dummy, NULL, FALSE, NULL);
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
    return betti[0];
    unsafe_wrap(Array, betti, (nrows, ncols), true);});


JULIA_CPP_MODULE_END
