#include "coeff_rings.h"

typedef struct __singular_coeff_ring_struct {
  int64_t has_simple_alloc;
  int64_t has_simple_Inverse;
  int64_t is_field;
  int64_t is_domain;
  int64_t ch;
  void *  data;
  void *  cfInit;
  void *  cfInt;
  void *  cfMPZ;
  void *  cfInpNeg;
  void *  cfCopy;
  void *  cfDelete;
  void *  cfAdd;
  void *  cfInpAdd;
  void *  cfSub;
  void *  cfMult;
  void *  cfInpMult;
  void *  cfDiv;
  void *  cfDivBy;
  void *  cfInvers;
  void *  cfAnn;
  void *  cfGcd;
  void *  cfSubringGcd;
  void *  cfExtGcd;
  void *  cfGreater;
  void *  cfEqual;
  void *  cfIsZero;
  void *  cfIsOne;
  void *  cfIsMOne;
  void *  cfGreaterZero;
  void *  cfWriteLong;
  void *  cfCoeffWrite;
  void *  cfFactorize;
} singular_coeff_ring_struct;

static BOOLEAN nCoeffIsEqual(const coeffs r, n_coeffType n, void *d)
{
  /* test, if r is an instance of nInitCoeffs(n,parameter) */
  return (n==r->type) &&(r->data==d);
}

void fill_coeffs_with_function_data(jl_value_t * coeff_struct, void * cf_void)
{

  coeffs cf = reinterpret_cast<coeffs>(cf_void);

  singular_coeff_ring_struct * cf_input =
      reinterpret_cast<singular_coeff_ring_struct *>(coeff_struct);
  cf->has_simple_Alloc = (BOOLEAN)cf_input->has_simple_alloc;
  cf->has_simple_Inverse = (BOOLEAN)cf_input->has_simple_Inverse;
  cf->is_field = (BOOLEAN)cf_input->is_field;
  cf->is_domain = (BOOLEAN)cf_input->is_domain;
  cf->ch = (int)cf_input->ch;
  cf->data = cf_input->data;
  cf->nCoeffIsEqual = nCoeffIsEqual;
  cf->cfInit = (number(*)(long, const coeffs))cf_input->cfInit;
  cf->cfInt = (long (*)(number &, const coeffs))cf_input->cfInt;
  cf->cfMPZ = (void (*)(__mpz_struct *, number &, const coeffs))cf_input->cfMPZ;
  cf->cfInpNeg = (number(*)(number, const coeffs))cf_input->cfInpNeg;
  cf->cfCopy = (number(*)(number, const coeffs))cf_input->cfCopy;
  cf->cfDelete = (void (*)(number *, const coeffs))cf_input->cfDelete;
  cf->cfAdd = (numberfunc)cf_input->cfAdd;
  cf->cfInpAdd = (void (*)(number &, number, const coeffs))cf_input->cfInpAdd;
  cf->cfSub = (numberfunc)cf_input->cfSub;
  cf->cfMult = (numberfunc)cf_input->cfMult;
  cf->cfInpMult = (void (*)(number &, number, const coeffs))cf_input->cfInpMult;
  cf->cfDiv = (numberfunc)cf_input->cfDiv;
  cf->cfDivBy = (BOOLEAN(*)(number, number, const coeffs))cf_input->cfDivBy;
  cf->cfInvers = (number(*)(number, const coeffs))cf_input->cfInvers;
  cf->cfAnn = (number(*)(number, const coeffs))cf_input->cfAnn;
  cf->cfGcd = (numberfunc)cf_input->cfGcd;
  cf->cfSubringGcd = (numberfunc)cf_input->cfSubringGcd;
  cf->cfExtGcd =
      (number(*)(number, number, number *, number *, const coeffs))cf_input->cfExtGcd;
  cf->cfGreater = (BOOLEAN(*)(number, number, const coeffs))cf_input->cfGreater;
  cf->cfEqual = (BOOLEAN(*)(number, number, const coeffs))cf_input->cfEqual;
  cf->cfIsZero = (BOOLEAN(*)(number, const coeffs))cf_input->cfIsZero;
  cf->cfIsOne = (BOOLEAN(*)(number, const coeffs))cf_input->cfIsOne;
  cf->cfIsMOne = (BOOLEAN(*)(number, const coeffs))cf_input->cfIsMOne;
  cf->cfGreaterZero = (BOOLEAN(*)(number, const coeffs))cf_input->cfGreaterZero;
  cf->cfWriteLong = (void (*)(number, const coeffs))cf_input->cfWriteLong;
  cf->cfCoeffWrite = (void (*)(const coeffs, BOOLEAN))cf_input->cfCoeffWrite;
  cf->cfFactorize = (nFactorizeFunc)cf_input->cfFactorize;
}

ideal copy_factorization_result(ideal factors, jlcxx::ArrayRef<int> exponents,
                                void * v_void, int with_exps,
                                ring source, void * target_void)
{
  ring target = reinterpret_cast<ring>(target_void);
  ideal result = idrCopyR(factors, source, target);

  if (with_exps != 1)
  {
    intvec ** v = reinterpret_cast<intvec **>(v_void);
    if (v != NULL)
    {
      const int n = si_max(1, IDELEMS(result));
      (*v) = new intvec(n);
      int * content = (*v)->ivGetVec();
      for (int i = 0; i < n; i++)
      {
        content[i] = (i < exponents.size()) ? exponents[i] : 1;
      }
    }
  }

  return result;
}


void singular_define_coeff_rings(jlcxx::Module & singular)
{
  singular.method("fill_coeffs_with_function_data", fill_coeffs_with_function_data);
  singular.method("nRegister", [](n_coeffType type, void * func) {
    return nRegister(type, (cfInitCharProc)func);
  });
  singular.method("get_coeff_data", [](coeffs c) {
    return c->data;
  });
  singular.method("get_coeff_data_void", [](void * c) {
    return reinterpret_cast<coeffs>(c)->data;
  });
  singular.method("cast_number_to_void", [](number n) {
    return reinterpret_cast<void *>(n);
  });
  singular.method("cast_void_to_number", [](void * n) {
    return reinterpret_cast<number>(n);
  });
  singular.method("factorization_callback_ring_copy", [](void * r) {
    return rCopy(reinterpret_cast<ring>(r));
  });
  singular.method("factorization_callback_poly_copy", [](void * p, void * r) {
    return p_Copy(reinterpret_cast<poly>(p), reinterpret_cast<ring>(r));
  });
  singular.method("factorization_callback_poly_copy_to_ring", [](void * p, void * source, ring target) {
    return prCopyR(reinterpret_cast<poly>(p), reinterpret_cast<ring>(source), target);
  });
  singular.method("factorization_callback_coeff_data", [](void * r) {
    return reinterpret_cast<ring>(r)->cf->data;
  });
  singular.method("copy_factorization_result", copy_factorization_result);
  singular.method("get_coeffs_BIGINT", []() {
    return coeffs_BIGINT;
  });
}
