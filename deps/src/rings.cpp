#include "rings.h"

auto rDefault_helper(coeffs                       cf,
                     jlcxx::ArrayRef<std::string> vars,
                     rRingOrder_t                 ord)
{
    auto    len = vars.size();
    char ** vars_ptr = new char *[len];
    for (int i = 0; i < len; i++) {
        vars_ptr[i] = new char[vars[i].length() + 1];
        std::strcpy(vars_ptr[i], vars[i].c_str());
    }
    auto r = rDefault(cf, len, vars_ptr, ord);
    delete[] vars_ptr;
    r->ShortOut = 0;
    return r;
}

auto rDefault_long_helper(coeffs                        cf,
                          jlcxx::ArrayRef<uint8_t *>    vars,
                          jlcxx::ArrayRef<rRingOrder_t> ord,
                          int *                         blk0,
                          int *                         blk1,
                          unsigned long                 bitmask)
{
    auto    len = vars.size();
    char ** vars_ptr = new char *[len];
    for (int i = 0; i < len; i++) {
        vars_ptr[i] = reinterpret_cast<char *>(vars[i]);
        // std::strcpy(vars_ptr[i],vars[i].c_str());
    }
    auto           len_ord = ord.size();
    rRingOrder_t * ord_ptr =
        (rRingOrder_t *)omAlloc0(len_ord * sizeof(rRingOrder_t));
    for (int i = 0; i < len_ord; i++) {
        ord_ptr[i] = ord[i];
    }
    int ** wvhdl = NULL;
    auto r = rDefault(cf, len, vars_ptr, len_ord, ord_ptr, blk0, blk1, wvhdl,
                      bitmask);
    delete[] vars_ptr;
    r->ShortOut = 0;
    return r;
}

void singular_define_rings(jlcxx::Module & Singular)
{
    Singular.method("rDefault_helper", &rDefault_helper);
    Singular.method("rDefault_long_helper", &rDefault_long_helper);
    Singular.method("rDelete", &rDelete);
    Singular.method("rString", [](ip_sring * r) {
        auto s = rString(r);
        return std::string(s);
    });
    Singular.method("rChar", &rChar);
    Singular.method("rGetVar", &rGetVar);
    Singular.method("rVar", &rVar);
    Singular.method("rGetExpSize", [](unsigned long bitmask, int N) {
        int bits;
        return static_cast<unsigned int>(rGetExpSize(bitmask, bits, N));
    });
    Singular.method("rHasGlobalOrdering", &rHasGlobalOrdering);
    Singular.method("rBitmask",
                    [](ip_sring * r) { return (unsigned int)r->bitmask; });
    Singular.method("p_Delete", [](spolyrec * p, ip_sring * r) {
        return p_Delete(&p, r);
    });
    Singular.method("p_Copy",
                    [](spolyrec * p, ip_sring * r) { return p_Copy(p, r); });
    Singular.method("p_IsOne",
                    [](spolyrec * p, ip_sring * r) { return p_IsOne(p, r); });
    Singular.method("p_One", [](ip_sring * r) { return p_One(r); });
    Singular.method("p_Init", [](ip_sring * r) { return p_Init(r); });
    Singular.method("p_IsUnit", [](spolyrec * p, ip_sring * r) {
        return p_IsUnit(p, r);
    });
    Singular.method("p_GetExp", [](spolyrec * p, int i, ip_sring * r) {
        return p_GetExp(p, i, r);
    });
    Singular.method("p_GetComp", [](spolyrec * p, ip_sring * r) {
        return p_GetComp(p, r);
    });
    Singular.method("p_String", [](spolyrec * p, ip_sring * r) {
        std::string s(p_String(p, r));
        return s;
    });
    Singular.method("p_ISet",
                    [](long i, ip_sring * r) { return p_ISet(i, r); });
    Singular.method("p_NSet",
                    [](snumber * p, ip_sring * r) { return p_NSet(p, r); });
    Singular.method("pLength", [](spolyrec * p) { return pLength(p); });
    Singular.method("SetpNext",
                    [](spolyrec * p, spolyrec * q) { p->next = q; });
    Singular.method("pNext", [](spolyrec * a) {
        poly p = pNext(a);
        return p;
    });
    Singular.method("p_Neg",
                    [](spolyrec * p, ip_sring * r) { return p_Neg(p, r); });
    Singular.method("pGetCoeff", [](spolyrec * p) { return pGetCoeff(p); });
    Singular.method("pSetCoeff", [](spolyrec * p, long c, ip_sring * r) {
        number n = n_Init(c, r);
        return p_SetCoeff(p, n, r);
    });
    Singular.method("pSetCoeff0", [](spolyrec * p, long c, ip_sring * r) {
        number n = n_Init(c, r);
        return p_SetCoeff0(p, n, r);
    });
    Singular.method("pLDeg", [](spolyrec * a, ip_sring * r) {
        long res;
        int  dummy;
        if (a != NULL) {
            res = r->pLDeg(a, &dummy, r);
        }
        else {
            res = -1;
        }
        return res;
    });
    Singular.method("p_Add_q", [](spolyrec * p, spolyrec * q, ip_sring * r) {
        return p_Add_q(p, q, r);
    });
    Singular.method("p_Sub", [](spolyrec * p, spolyrec * q, ip_sring * r) {
        return p_Sub(p, q, r);
    });
    Singular.method("p_Mult_q", [](spolyrec * p, spolyrec * q, ip_sring * r) {
        return p_Mult_q(p, q, r);
    });
    Singular.method("p_Power", [](spolyrec * p, int q, ip_sring * r) {
        return p_Power(p, q, r);
    });
    Singular.method("p_EqualPolys",
                    [](spolyrec * p, spolyrec * q, ip_sring * r) {
                        return p_EqualPolys(p, q, r);
                    });
    Singular.method("p_Divide", [](spolyrec * p, spolyrec * q, ip_sring * r) {
        return p_Divide(p, q, r);
    });
    Singular.method("singclap_gcd",
                    [](spolyrec * p, spolyrec * q, ip_sring * r) {
                        return singclap_gcd(p, q, r);
                    });
    Singular.method("p_ExtGcd_internal", [](spolyrec * a, spolyrec * b,
                                            void * res, void * s, void * t,
                                            ip_sring * r) {
        return singclap_extgcd(a, b, reinterpret_cast<spolyrec *&>(res),
                               reinterpret_cast<spolyrec *&>(s),
                               reinterpret_cast<spolyrec *&>(t), r);
    });
    Singular.method("p_Content", [](spolyrec * p, ip_sring * r) {
        return p_Content(p, r);
    });
    Singular.method("p_GetExpVL_internal",
                    [](spolyrec * p, long * ev, ip_sring * r) {
                        return p_GetExpVL(p, ev, r);
                    });
    Singular.method("p_SetExpV_internal",
                    [](spolyrec * p, int * ev, ip_sring * r) {
                        return p_SetExpV(p, ev, r);
                    });
    Singular.method("p_Reduce",
                    [](spolyrec * p, sip_sideal * G, ip_sring * R) {
                        const ring origin = currRing;
                        rChangeCurrRing(R);
                        poly res = kNF(G, R->qideal, p);
                        rChangeCurrRing(origin);
                        return res;
                    });
    Singular.method("p_Reduce",
                    [](sip_sideal * p, sip_sideal * G, ip_sring * R) {
                        const ring origin = currRing;
                        rChangeCurrRing(R);
                        ideal res = kNF(G, R->qideal, p);
                        rChangeCurrRing(origin);
                        return res;
                    });

    Singular.method("letterplace_ring_helper",
                    [](ip_sring * r, long block_size) {
                        rUnComplete(r);
                        r->isLPring = block_size;
                        r->ShortOut = FALSE;
                        r->CanShortOut = FALSE;
                        rComplete(r);
                    });
}
