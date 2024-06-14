#include "coeffs.h"

auto transExt_helper(coeffs cf, jlcxx::ArrayRef<uint8_t *> param)
{
  auto    len = param.size();
  char ** param_ptr = new char *[len];
  for (int i = 0; i < len; i++)
  {
    param_ptr[i] = reinterpret_cast<char *>(param[i]);
  }
  // use lex because the ordering doesn't matter for coefficient operations
  ring r = rDefault(cf, len, param_ptr, ringorder_lp);
  delete[] param_ptr;
  TransExtInfo extParam;
  extParam.r = r;
  return nInitChar(n_transExt, &extParam);
}

poly transExt_to_poly(number a, coeffs cf, ring r)
{
  // zero is represented by nullptr
  if (a == NULL || NUM((fraction)a) == NULL)
    return NULL;

  assume(cf->extRing != NULL);
  ring       ext = cf->extRing;
  nMapFunc   nMap = n_SetMap(ext->cf, r->cf);
  const ring origin = currRing;
  rChangeCurrRing(r);
  poly p = p_PermPoly(NUM((fraction)a), NULL, ext, r, nMap, NULL);
  rChangeCurrRing(origin);
  return p;
}

void singular_define_coeffs(jlcxx::Module & Singular)
{
  /* initialise a coefficient ring */
  Singular.method("nInitChar", &nInitChar);
  /* Helper to construct transcendental Extensions */
  Singular.method("transExt_helper", &transExt_helper);
  /* Helper to convert n_transExt into a polynomial */
  Singular.method("transExt_to_poly", &transExt_to_poly);
  /* get the characteristic of a coefficient domain */
  Singular.method("n_GetChar", [](coeffs n) {
    return n_GetChar(n);
  });

  /*
      given F = R(x), a univariate function field over R = QQ or Fp,
          construct the number field R[x]/minpoly.
      Return a new coefficient ring or a copy of the old one in case of error.
      Thus, the return should be checked to see if it is a number field
          (nCoeff_is_algExt) before use.

      TODO: this will be in the kernel soon
  */
  Singular.method("transExt_SetMinpoly", [](coeffs F, number a) {
    if (!nCoeff_is_transExt(F) || rVar(F->extRing) != 1)
    {
      WerrorS("cannot set minpoly for these coeffients");
      return nCopyCoeff(F);
    }

    BOOLEAN redefine_from_algext = FALSE;

    number p = n_Copy(a, F);
    n_Normalize(p, F);

    if (n_IsZero(p, F))
    {
      n_Delete(&p, F);
      return nCopyCoeff(F);
    }

    AlgExtInfo A;

    A.r = rCopy(F->extRing);    // Copy  ground field!
    // if minpoly was already set:
    if (F->extRing->qideal != NULL)
      id_Delete(&A.r->qideal, A.r);
    ideal q = idInit(1, 1);
    if (p == NULL || NUM((fraction)p) == NULL)
    {
      WerrorS("Could not construct the alg. extension: minpoly==0");
      // cleanup A: TODO
      rDelete(A.r);
      return nCopyCoeff(F);
    }
    if (!redefine_from_algext &&
        (DEN((fraction)(p)) !=
         NULL))    // minpoly must be a fraction with poly numerator...!!
    {
      poly n = DEN((fraction)(p));
      if (!p_IsConstantPoly(n, F->extRing))
      {
        WarnS("denominator must be constant - ignoring it");
      }
      p_Delete(&n, F->extRing);
      DEN((fraction)(p)) = NULL;
    }

    if (redefine_from_algext)
      q->m[0] = (poly)p;
    else
      q->m[0] = NUM((fraction)p);
    A.r->qideal = q;

    // :(
    //  NUM((fractionObject *)p) = NULL; // makes 0/ NULL fraction - which should not
    //  happen! n_Delete(&p, r->F); // doesn't expect 0/ NULL :(
    if (!redefine_from_algext)
    {
      EXTERN_VAR omBin fractionObjectBin;
      NUM((fractionObject *)p) = NULL;    // not necessary, but still...
      omFreeBin((ADDRESS)p, fractionObjectBin);
    }

    coeffs K = nInitChar(n_algExt, &A);
    if (K == NULL)
    {
      WerrorS("Could not construct the alg. extension: llegal minpoly?");
      // cleanup A: TODO
      rDelete(A.r);
      return nCopyCoeff(F);
    }
    else
    {
      return K;
    }
  });

  /*
      if K = R[x]/minpoly(x) (with nCoeff_is_algExt(cf) = true),
      return the corresponding transcendental extension R(x): this is the
      parent of minpoly.
  */
  Singular.method("algExt_UnsetMinpoly", [](coeffs K) {
    if (nCoeff_is_algExt(K) && !nCoeff_is_GF(K))
    {
      ring K0 = rCopy0(K->extRing, FALSE, TRUE);
      rComplete(K0);
      TransExtInfo e;
      e.r = K0;
      coeffs F = nInitChar(n_transExt, &e);
      return F;
    }
    else
    {
      WerrorS("cannot unset minpoly for these coeffients");
      return nCopyCoeff(K);
    }
  });

  /*
      Given K = R[x]/minpoly(x) and the corresponding F = R(x) return
      the minpoly as an element of F.
  */
  Singular.method("algExt_GetMinpoly", [](coeffs K, coeffs F) {
    if (nCoeff_is_algExt(K) && !nCoeff_is_GF(K))
    {
      nMapFunc nMap = n_SetMap(K, F);
      number   minpoly = reinterpret_cast<number>(K->extRing->qideal->m[0]);
      return nMap(minpoly, K, F);
    }
    else
    {
      WerrorS("cannot get minpoly for these coeffients");
      return n_Init(0, F);
    }
  });

  /*
      Given an element a of K = R[x]/minpoly(x) and the corresponding
      F = R(x), return the element a as an element of F.
  */
  Singular.method("algExt_to_transExt", [](number a, coeffs K, coeffs F) {
    if (nCoeff_is_algExt(K) && !nCoeff_is_GF(K))
    {
      nMapFunc nMap = n_SetMap(K, F);
      return nMap(a, K, F);
    }
    else
    {
      WerrorS("cannot use algExt_to_transExt for these coeffients");
      return n_Init(0, F);
    }
  });

  /*
      Given an element a of F = R(x) and K = R[x]/minpoly(x),
      return a as an element of K.
  */
  Singular.method("transExt_to_algExt", [](number a, coeffs K, coeffs F) {
    if (nCoeff_is_algExt(K) && !nCoeff_is_GF(K))
    {
      // naCopyTrans2AlgExt has a bug with zero input
      if (a == NULL || NUM((fraction)a) == NULL)
        return reinterpret_cast<number>(NULL);

      nMapFunc nMap = n_SetMap(F, K);
      return nMap(a, F, K);
    }
    else
    {
      WerrorS("cannot use transExt_to_algExt for these coeffients");
      return n_Init(0, F);
    }
  });

  /* only used to get the order of a N_GField */
  Singular.method("nfCharQ", [](coeffs n) {
    return n->m_nfCharQ;
  });

  Singular.method("nCoeff_is_Zp", [](coeffs n) {
    return bool(nCoeff_is_Zp(n));
  });

  Singular.method("nCoeff_is_Q", [](coeffs n) {
    return bool(nCoeff_is_Q(n));
  });

  Singular.method("nCoeff_is_Z", [](coeffs n) {
    return bool(nCoeff_is_Z(n));
  });

  Singular.method("nCoeff_is_GF", [](coeffs n) {
    return bool(nCoeff_is_GF(n));
  });

  Singular.method("nCoeff_is_transExt", [](coeffs n) {
    return bool(nCoeff_is_transExt(n));
  });

  Singular.method("nCoeff_is_algExt", [](coeffs n) {
    return bool(nCoeff_is_algExt(n));
  });

  Singular.method("nCoeff_is_Nemo_AnticNumberField", [](coeffs n) {
    return n->type==n_Nemo_AnticNumberField;
  });

  Singular.method("nCoeff_is_Nemo_QQField", [](coeffs n) {
    return n->type==n_Nemo_QQField;
  });

  Singular.method("nCoeff_is_Nemo_ZZRing", [](coeffs n) {
    return n->type==n_Nemo_ZZRing;
  });

  Singular.method("nCoeff_is_Nemo_FqPolyRepField", [](coeffs n) {
    return n->type==n_Nemo_FqPolyRepField;
  });

  Singular.method("nCoeff_is_Nemo_fqPolyRepField", [](coeffs n) {
    return n->type==n_Nemo_fqPolyRepField;
  });

  Singular.method("nCoeff_is_Nemo_Field", [](coeffs n) {
    return n->type==n_Nemo_Field;
  });

  Singular.method("nCoeff_is_Nemo_Ring", [](coeffs n) {
    return n->type==n_Nemo_Ring;
  });

  Singular.method("nGetCoeffData", [](coeffs n) {
    return n->data;
  });

  /* make a copy of a coefficient domain (actually just increments a
   * reference count) */
  Singular.method("nCopyCoeff", &nCopyCoeff);
  /* kill a coefficient ring */
  Singular.method("nKillChar", &nKillChar);
  /* return a function to convert between rings */
  Singular.method("n_SetMap", [](const coeffs x, const coeffs y) {
    return reinterpret_cast<void *>(n_SetMap(x, y));
  });

  /* Identity function on coefficient ring */
  Singular.method("ndCopyMap", []() {
    return reinterpret_cast<void *>(ndCopyMap);
  });

  Singular.method("nApplyMapFunc", [](void * map, snumber * x, coeffs a, coeffs b) {
    return reinterpret_cast<nMapFunc>(map)(x, a, b);
  });

  Singular.method("n_Init", [](long x, coeffs n) {
    return n_Init(x, n);
  });

  Singular.method("n_Copy", [](snumber * x, const coeffs n) {
    return n_Copy(x, n);
  });

  Singular.method("nCoeff_has_simple_Alloc", [](coeffs x) {
    return nCoeff_has_simple_Alloc(x) > 0;
  });

  Singular.method("n_GetMPZ_internal", [](void * ptr, number n, coeffs x) {
    n_MPZ(reinterpret_cast<__mpz_struct *>(ptr), n, x);
  });

  Singular.method("n_InitMPZ_internal", [](void * ptr, coeffs x) {
    return n_InitMPZ(reinterpret_cast<__mpz_struct *>(ptr), x);
  });

  Singular.method("n_Delete", [](snumber * n, coeffs cf) {
    if (n != NULL)
    {
      n_Delete(&n, cf);
    }
  });

  Singular.method("n_Write_internal", [](snumber * x, coeffs cf, const int d) {
    return n_Write(x, cf, d);
  });

  Singular.method("n_Add", [](snumber * a, snumber * b, coeffs c) {
    return n_Add(a, b, c);
  });

  Singular.method("n_Sub", [](snumber * a, snumber * b, coeffs c) {
    return n_Sub(a, b, c);
  });

  Singular.method("n_Mult", [](snumber * a, snumber * b, coeffs c) {
    return n_Mult(a, b, c);
  });

  Singular.method("n_Neg", [](snumber * a, coeffs c) {
    number nn = n_Copy(a, c);
    nn = n_InpNeg(nn, c);
    return nn;
  });

  Singular.method("n_Invers", [](snumber * a, coeffs c) {
    return n_Invers(a, c);
  });

  Singular.method("n_ExactDiv", [](snumber * a, snumber * b, coeffs c) {
    return n_ExactDiv(a, b, c);
  });

  Singular.method("n_Div", [](snumber * a, snumber * b, coeffs c) {
    number z = n_Div(a, b, c);
    n_Normalize(z, c);
    return z;
  });

  /* return the numerator of a, possibly modifying a in the process */
  Singular.method("n_GetNumerator", [](number & a, coeffs c) {
    return n_GetNumerator(a, c);
  });

  /* return the denominator of a, possibly modifying a in the process */
  Singular.method("n_GetDenom", [](number & a, coeffs c) {
    return n_GetDenom(a, c);
  });

  /*
      since references are such a pain on the julia side,
      n_Normalize takes a non-reference and returns the modifed result
  */
  Singular.method("n_Normalize", [](number a, coeffs c) {
    n_Normalize(a, c);
    return a;
  });

  Singular.method("n_Power", [](snumber * a, int b, coeffs c) {
    number res;
    n_Power(a, b, &res, c);
    return res;
  });

  Singular.method("n_Gcd", [](snumber * a, snumber * b, coeffs c) {
    return n_Gcd(a, b, c);
  });

  Singular.method("n_SubringGcd", [](snumber * a, snumber * b, coeffs c) {
    return n_SubringGcd(a, b, c);
  });

  Singular.method("n_Lcm", [](snumber * a, snumber * b, coeffs c) {
    return n_Lcm(a, b, c);
  });

  Singular.method("n_ExtGcd",
                  [](number a, number b, number * s, number * t, const coeffs c) {
    return n_ExtGcd(a, b, s, t, c);
  });

  Singular.method("n_IsZero", [](snumber * x, const coeffs n) {
    return n_IsZero(x, n) > 0;
  });

  Singular.method("n_IsOne", [](snumber * x, const coeffs n) {
    return n_IsOne(x, n) > 0;
  });

  Singular.method("n_Greater", [](snumber * x, snumber * y, const coeffs n) {
    return n_Greater(x, y, n) > 0;
  });

  Singular.method("n_GreaterZero", [](snumber * x, const coeffs n) {
    return n_GreaterZero(x, n) > 0;
  });

  Singular.method("n_Equal", [](snumber * x, snumber * y, const coeffs n) {
    return n_Equal(x, y, n) > 0;
  });

  Singular.method("n_InpAdd", [](snumber * x, snumber * y, const coeffs n) {
    n_InpAdd(x, y, n);
    return x;
  });

  Singular.method("n_InpMult", [](snumber * x, snumber * y, const coeffs n) {
    n_InpMult(x, y, n);
    return x;
  });

  Singular.method("n_QuotRem", [](number x, number y, number * r, const coeffs n) {
    return n_QuotRem(x, y, r, n);
  });

  Singular.method("n_IntMod", &n_IntMod);

  Singular.method("n_Farey", &n_Farey);

  Singular.method("n_ChineseRemainderSym",
                  [](void * x, void * y, int n, int sign_flag, coeffs c) {
    CFArray inv_cache(n);
    return n_ChineseRemainderSym(reinterpret_cast<snumber **>(x),
                                 reinterpret_cast<snumber **>(y), n, sign_flag, inv_cache,
                                 c);
  });

  Singular.method("n_Param", [](int x, const coeffs n) {
    return n_Param(x, n);
  });

  Singular.method("n_NumberOfParameters", [](coeffs r) {
    return n_NumberOfParameters(r);
  });

  Singular.method("n_ParameterName", [](int i, coeffs r) {
    return std::string(n_ParameterNames(r)[i]);
  });

  Singular.method("StringSetS_internal", [](std::string m) {
    return StringSetS(m.c_str());
  });

  Singular.method("StringEndS", []() {
    char *      m = StringEndS();
    std::string s(m);
    omFree(m);
    return s;
  });

  Singular.method("omAlloc0", [](size_t siz) {
    return (void *)omAlloc0(siz);
  });

  Singular.method("omFree_internal", [](void * m) {
    omFree(m);
  });

  /* Setting a Ptr{number} to a number */
  Singular.method("setindex_internal", [](void * x, snumber * y) {
    *reinterpret_cast<snumber **>(x) = y;
  });

  Singular.method("setindex_internal_void", [](void * x, void * y) {
    reinterpret_cast<void **>(x)[0] = y;
  });

  Singular.method("mpz_init_set_internal", [](void * x, void * y) {
    mpz_init_set(reinterpret_cast<mpz_ptr>(x), reinterpret_cast<mpz_ptr>(y));
  });

  Singular.method("mpz_init_set_si_internal", [](void * x, long y) {
    mpz_init_set_si(reinterpret_cast<mpz_ptr>(x), y);
  });
}
