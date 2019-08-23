#include "ideals.h"

auto id_sres_helper(sip_sideal * m, int n, ring R)
{
    auto origin = currRing;
    rChangeCurrRing(R);
    syStrategy s = sySchreyer(m, n);
    rChangeCurrRing(origin);
    auto r = s->minres;
    bool minimal = true;
    if (r == NULL) {
        r = s->fullres;
        minimal = false;
    }
    return std::make_tuple(reinterpret_cast<void *>(r), s->length, minimal);
}


auto id_fres_helper(sip_sideal * I, int n, std::string method, ring R)
{
    auto origin = currRing;
    rChangeCurrRing(R);
    syStrategy s = syFrank(I, n, method.c_str());
    rChangeCurrRing(origin);
    auto r = s->minres;
    bool minimal = true;
    if (r == NULL) {
        r = s->fullres;
        minimal = false;
    }
    return std::make_tuple(reinterpret_cast<void *>(r), s->length, minimal);
}


ideal id_Syzygies_internal(ideal m, ring o)
{
    ideal      id = NULL;
    intvec *   n = NULL;
    tHomog     h = testHomog;
    const ring origin = currRing;
    rChangeCurrRing(o);
    id = idSyzygies(m, h, &n);
    rChangeCurrRing(origin);
    if (n != NULL)
        delete n;
    return id;
}

auto id_Slimgb_helper(ideal a, ring b, bool complete_reduction = false)
{
    //  bool complete_reduction= false;
    unsigned int crbit;
    if (complete_reduction)
        auto crbit = Sy_bit(OPT_REDSB);
    else
        crbit = 0;
    ideal id = NULL;
    if (!idIs0(a)) {
        intvec *     n = NULL;
        tHomog       h = testHomog;
        const ring   origin = currRing;
        unsigned int save_opt = si_opt_1;
        si_opt_1 |= crbit;
        rChangeCurrRing(b);
        id = t_rep_gb(b, a, a->rank);
        si_opt_1 = save_opt;
        rChangeCurrRing(origin);
        if (n != NULL)
            delete n;
    }
    else
        id = idInit(0, a->rank);
    return id;
}

auto id_Std_helper(ideal a, ring b, bool complete_reduction = false)
{
    // bool complete_reduction= false;
    unsigned int crbit;
    if (complete_reduction)
        crbit = Sy_bit(OPT_REDSB);
    else
        crbit = 0;
    ideal id = NULL;
    if (!idIs0(a)) {
        intvec *     n = NULL;
        tHomog       h = testHomog;
        const ring   origin = currRing;
        unsigned int save_opt = si_opt_1;
        si_opt_1 |= crbit;
        rChangeCurrRing(b);
        id = kStd(a, b->qideal, h, &n);
        si_opt_1 = save_opt;
        rChangeCurrRing(origin);
        if (n != NULL)
            delete n;
    }
    else
        id = idInit(0, a->rank);
    return id;
}

void singular_define_ideals(jlcxx::Module & Singular)
{
    Singular.method("id_Delete",
                    [](ideal m, ring n) { return id_Delete(&m, n); });

    Singular.method("id_Copy", &id_Copy);

    Singular.method("idInit", &idInit);

    Singular.method("setindex_internal",
                    [](ideal r, poly n, int o) { return r->m[o] = n; });

    Singular.method("getindex",
                    [](ideal r, int o) { return (poly)(r->m[o]); });

    Singular.method("idIs0", &idIs0);

    Singular.method("id_IsConstant", &id_IsConstant);

    Singular.method("id_IsZeroDim", &id_IsZeroDim);

    Singular.method("idElem", &idElem);

    Singular.method("id_Normalize", &id_Normalize);

    Singular.method("id_Head", &id_Head);

    Singular.method("id_MaxIdeal",
                    [](int m, ring n) { return id_MaxIdeal(m, n); });

    Singular.method("id_Add", &id_Add);

    Singular.method("id_Mult", &id_Mult);

    Singular.method("id_Power", &id_Power);

    Singular.method("id_IsEqual", [](ideal m, ideal n, ring o) {
        return mp_Equal((ip_smatrix *)m, (ip_smatrix *)n, o);
    });

    Singular.method("id_FreeModule", &id_FreeModule);

    Singular.method("idSkipZeroes", &idSkipZeroes);

    Singular.method("ngens", [](ideal m) { return (int)IDELEMS(m); });

    Singular.method("rank", [](ideal m) { return (int)m->rank; });

    Singular.method("id_Quotient", [](ideal a, ideal b, bool c, ring d) {
        const ring origin = currRing;
        rChangeCurrRing(d);
        ideal id = idQuot(a, b, c, TRUE);
        rChangeCurrRing(origin);
        return id;
    });

    Singular.method("id_Intersection", [](ideal a, ideal b, ring c) {
        const ring origin = currRing;
        rChangeCurrRing(c);
        ideal id = idSect(a, b);
        rChangeCurrRing(origin);
        return id;
    });

    Singular.method("id_Syzygies", &id_Syzygies_internal);

    Singular.method("id_sres", &id_sres_helper);

    Singular.method("id_fres", &id_fres_helper);

    Singular.method("id_Slimgb", &id_Slimgb_helper);

    Singular.method("id_Std", &id_Std_helper);

    Singular.method("id_Eliminate", [](ideal m, poly p, ring o) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal res = idElimination(m, p);
        rChangeCurrRing(origin);
        return res;
    });

    Singular.method("id_Satstd", &id_Satstd);

    Singular.method("id_Array2Vector", [](void * p, int a, ring o) {
        return id_Array2Vector(reinterpret_cast<poly *>(p), a, o);
    });

    Singular.method("p_Vector2Array", [](poly p, void * a, int b, ring o) {
        p_Vec2Array(p, reinterpret_cast<poly *>(a), b, o);
    });

    Singular.method("internal_void_to_poly_helper",
                    [](void * p) { return reinterpret_cast<poly>(p); });

    Singular.method(
        "maGetPreimage", [](ring trgt, ideal a, ideal b, ring src) {
            sip_smap sing_map = {a->m, (char *)"julia_ring", 1, a->ncols};
            return maGetPreimage(trgt, &sing_map, b, src);
    });

    Singular.method("id_Jet", [](ideal I, int n, ring r) {
        ideal res = id_Jet(I, n, r);
        return res;
    });

    Singular.method("id_vdim", [](ideal I, ring r) {
        const ring origin = currRing;
        rChangeCurrRing(r);
	int n=scMult0Int(I, r->qideal);
	rChangeCurrRing(origin);
	return n;
    });

    Singular.method("id_kbase", [](ideal I, ring r) {
        ideal res;
	const ring origin = currRing;
        rChangeCurrRing(r);
        res = scKBase(-1, I, r->qideal);
        rChangeCurrRing(origin);
        return res;
    });

    Singular.method("id_highcorner", [](ideal I, ring r) {
        poly h;
        const ring origin = currRing;
        rChangeCurrRing(r);
        h = iiHighCorner(I, 0);
        rChangeCurrRing(origin);
        return h;
    });
    Singular.method("maMapIdeal", [](ideal map_id, ring pr, ideal im_id,
                    ring im) {

        rChangeCurrRing(pr);
        nMapFunc nMap =n_SetMap(currRing->cf, im->cf);
        return maMapIdeal(map_id, pr, im_id, im, nMap);
    });
    Singular.method("idMinBase", [](ideal I, ring r) {
        rChangeCurrRing(r);
        return idMinBase(I);
    });
}
