#include "coeffs.h"

auto transExt_helper(coeffs cf, jlcxx::ArrayRef<uint8_t *> param)
{
    auto    len = param.size();
    char ** param_ptr = new char *[len];
    for (int i = 0; i < len; i++) {
        param_ptr[i] = reinterpret_cast<char *>(param[i]);
    }
    ring r = rDefault(cf, len, param_ptr);
    r->order[0] = ringorder_dp;
    delete[] param_ptr;
    TransExtInfo extParam;
    extParam.r = r;
    return nInitChar(n_transExt, &extParam);
}

auto transExt_to_poly(number a, coeffs cf, ring r)
{
    assume(cf->extRing != NULL);
    ring ext = cf->extRing;
    fraction f = (fraction)a;
    nMapFunc nMap = n_SetMap(ext->cf, r->cf);
    rChangeCurrRing(r);
    return p_PermPoly(f->numerator, NULL, ext, r, nMap, NULL);
}
void singular_define_coeffs(jlcxx::Module & Singular)
{
    /* initialise a coefficient ring */
    Singular.method("nInitChar", &nInitChar);
    /* Helper to construct transcendental Extensions */
    Singular.method("transExt_helper", &transExt_helper);
    Singular.method("transExt_to_poly", &transExt_to_poly);
    /* get the characteristic of a coefficient domain */
    Singular.method("n_GetChar", [](coeffs n) { return n_GetChar(n); });

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

    Singular.method("nApplyMapFunc",
                    [](void * map, snumber * x, coeffs a, coeffs b) {
                        return reinterpret_cast<nMapFunc>(map)(x, a, b);
                    });

    Singular.method("n_Init", [](long x, coeffs n) { return n_Init(x, n); });

    Singular.method("n_Copy",
                    [](snumber * x, const coeffs n) { return n_Copy(x, n); });

    Singular.method("nCoeff_has_simple_Alloc",
                    [](coeffs x) { return nCoeff_has_simple_Alloc(x) > 0; });

    Singular.method("n_GetMPZ_internal", [](void * ptr, number n, coeffs x) {
        n_MPZ(reinterpret_cast<__mpz_struct *>(ptr), n, x);
    });

    Singular.method("n_InitMPZ_internal", [](void * ptr, coeffs x) {
        return n_InitMPZ(reinterpret_cast<__mpz_struct *>(ptr), x);
    });

    Singular.method("n_Delete", [](snumber * n, coeffs cf) {
        if (n != NULL) {
            n_Delete(&n, cf);
        }
    });

    Singular.method("n_Delete_Q", [](snumber * n, coeffs cf) {
        if (n != NULL) {
            n_Delete(&n, cf);
        }
    });

    Singular.method("n_Write_internal",
                    [](snumber * x, coeffs cf, const int d) {
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

    Singular.method("n_Invers",
                    [](snumber * a, coeffs c) { return n_Invers(a, c); });

    Singular.method("n_ExactDiv", [](snumber * a, snumber * b, coeffs c) {
        return n_ExactDiv(a, b, c);
    });

    Singular.method("n_Div", [](snumber * a, snumber * b, coeffs c) {
        number z = n_Div(a, b, c);
        n_Normalize(z, c);
        return z;
    });

    Singular.method("n_GetNumerator", [](snumber * a, coeffs c) {
        return n_GetNumerator(a, c);
    });

    Singular.method("n_GetDenom",
                    [](snumber * a, coeffs c) { return n_GetDenom(a, c); });

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

    Singular.method("n_ExtGcd_internal", [](snumber * a, snumber * b,
                                            void * s, void * t, coeffs c) {
        return n_ExtGcd(a, b, reinterpret_cast<snumber **>(s),
                        reinterpret_cast<snumber **>(t), c);
    });

    Singular.method("n_IsZero", [](snumber * x, const coeffs n) {
        return n_IsZero(x, n) > 0;
    });

    Singular.method("n_IsOne", [](snumber * x, const coeffs n) {
        return n_IsOne(x, n) > 0;
    });

    Singular.method("n_Greater",
                    [](snumber * x, snumber * y, const coeffs n) {
                        return n_Greater(x, y, n) > 0;
                    });

    Singular.method("n_GreaterZero", [](snumber * x, const coeffs n) {
        return n_GreaterZero(x, n) > 0;
    });

    Singular.method("n_Equal", [](snumber * x, snumber * y, const coeffs n) {
        return n_Equal(x, y, n) > 0;
    });

    Singular.method("n_InpAdd", [](snumber * x, snumber * y, const coeffs n) {
        snumber * xx = x;
        n_InpAdd(xx, y, n);
        return xx;
    });

    Singular.method("n_InpMult",
                    [](snumber * x, snumber * y, const coeffs n) {
                        return n_InpMult(x, y, n);
                    });

    Singular.method("n_QuotRem_internal", [](snumber * x, snumber * y,
                                             void * p, const coeffs n) {
        return n_QuotRem(x, y, reinterpret_cast<snumber **>(p), n);
    });

    Singular.method("n_Rem", [](snumber * x, snumber * y, const coeffs n) {
        number qq;
        return n_QuotRem(x, y, &qq, n);
    });

    Singular.method("n_IntMod", &n_IntMod);

    Singular.method("n_Farey", &n_Farey);

    Singular.method("n_ChineseRemainderSym", [](void * x, void * y, int n,
                                                int sign_flag, coeffs c) {
        CFArray inv_cache(n);
        return n_ChineseRemainderSym(reinterpret_cast<snumber **>(x),
                                     reinterpret_cast<snumber **>(y), n,
                                     sign_flag, inv_cache, c);
    });

    Singular.method("n_Param",
                    [](int x, const coeffs n) { return n_Param(x, n); });

    Singular.method("StringSetS_internal",
                    [](std::string m) { return StringSetS(m.c_str()); });

    Singular.method("StringEndS", []() {
        char *      m = StringEndS();
        std::string s(m);
        omFree(m);
        return s;
    });

    Singular.method("omAlloc0",
                    [](size_t siz) { return (void *)omAlloc0(siz); });

    Singular.method("omFree_internal", [](void * m) { omFree(m); });

    /* Setting a Ptr{number} to a number */

    Singular.method("setindex_internal", [](void * x, snumber * y) {
        *reinterpret_cast<snumber **>(x) = y;
    });

    Singular.method("setindex_internal_void", [](void * x, void * y) {
        reinterpret_cast<void **>(x)[0] = y;
    });

    Singular.method("mpz_init_set_internal", [](void* x, void* y ){
        mpz_init_set(reinterpret_cast<mpz_ptr>(x),reinterpret_cast<mpz_ptr>(y));
    });

    Singular.method("mpz_init_set_si_internal", [](void* x, long y ){
        mpz_init_set_si(reinterpret_cast<mpz_ptr>(x),y);
    });
}
