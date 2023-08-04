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
    return std::make_tuple(s, minimal);
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
    return std::make_tuple(s, minimal);
}

auto id_res_helper(sip_sideal * I, int n, int minimize, ring R)
{
    auto origin = currRing;
    rChangeCurrRing(R);
    syStrategy s = syResolution(I, n, NULL, (BOOLEAN)minimize);
    rChangeCurrRing(origin);
    auto r = s->minres;
    bool minimal = true;
    if (r == NULL) {
        r = s->fullres;
        minimal = false;
    }
    for(int i=0;i<=n+1;i++)
    {
      if (r[i]==NULL)
      {
        r[i]=idInit(1,1);
        break;
      }
    }
    return std::make_tuple(s, minimal);
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

auto id_InterRed_helper(ideal a, ring b)
{
    ideal id = NULL;
    if (!idIs0(a)) {
        const ring   origin = currRing;
        rChangeCurrRing(b);
        id = kInterRed(a, b->qideal);
        rChangeCurrRing(origin);
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

auto id_MinStd_helper(ideal a, ring b, bool complete_reduction = false)
{
    // bool complete_reduction= false;
    unsigned int crbit;
    if (complete_reduction)
        crbit = Sy_bit(OPT_REDSB);
    else
        crbit = 0;
    ideal id = NULL;
    ideal m = NULL;
    if (!idIs0(a))
    {
        tHomog       h = testHomog;
        const ring   origin = currRing;
        unsigned int save_opt = si_opt_1;
        si_opt_1 |= crbit;
        rChangeCurrRing(b);
        id = kMin_std(a, b->qideal, h, NULL, m);
        si_opt_1 = save_opt;
        rChangeCurrRing(origin);
    }
    else
    {
        id = idInit(0, a->rank);
        m = idInit(0, a->rank);
    }
    return std::make_tuple(id, m);
}

intvec* to_intvec(jlcxx::ArrayRef<int> a)
{
    int sz = a.size();
    intvec * w = new intvec(sz);
    int * hi = w->ivGetVec();
    for (int i=0; i<sz; i++)
        hi[i] = a[i];
    return w;
}

auto id_StdHilb_helper(ideal a, ring b, jlcxx::ArrayRef<int> h, bool complete_reduction = false)
{
    intvec* hilb = to_intvec(h);

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
        id = kStd(a, b->qideal, h, &n, hilb);
        si_opt_1 = save_opt;
        rChangeCurrRing(origin);
        if (n != NULL)
            delete n;
    }
    else
        id = idInit(0, a->rank);
    delete hilb;
    return id;
}

auto id_StdHilbWeighted_helper(ideal a, ring b, jlcxx::ArrayRef<int> h, jlcxx::ArrayRef<int> vw, bool complete_reduction = false)
{
    intvec* hilb = to_intvec(h);
    intvec* varweights = to_intvec(vw);

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
        id = kStd(a,
              currRing->qideal,
              h,
              &n,           // module weights
              hilb,         // hilbert series
              0,0,          // syzComp, newIdeal
              varweights);  // weights of vars
        si_opt_1 = save_opt;
        rChangeCurrRing(origin);
        if (n != NULL)
            delete n;
    }
    else
        id = idInit(0, a->rank);
    delete hilb;
    delete varweights;
    return id;
}

auto id_TwoStd_helper(ideal a, ring b)
{
    const ring origin = currRing;
    rChangeCurrRing(b);
    ideal id = twostd(a);
    rChangeCurrRing(origin);
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

    Singular.set_override_module(jl_base_module);
    Singular.method("getindex",
                    [](ideal r, int o) { return (poly)(r->m[o]); });
    Singular.unset_override_module();

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

    Singular.method("id_MultP", [](ideal i, poly p, ring r) {
        return (ideal) mp_MultP((matrix) i, p, r);
    });

    Singular.method("pMultId", [](poly p, ideal i, ring r) {
        return (ideal) pMultMp(p, (matrix) i, r);
    });

    Singular.method("id_Power", &id_Power);

    Singular.method("id_IsEqual", [](ideal m, ideal n, ring o) {
        return mp_Equal((ip_smatrix *)m, (ip_smatrix *)n, o);
    });

    Singular.method("id_FreeModule", &id_FreeModule);

    Singular.method("idSkipZeroes", &idSkipZeroes);

    Singular.method("ngens", [](ideal m) { return (int)IDELEMS(m); });

    Singular.method("rank", [](ideal m) { return (int)m->rank; });

    Singular.method("id_Homogenize", &id_Homogenize);

    Singular.method("id_HomogenizeW", [](ideal a, int v, jlcxx::ArrayRef<int> w, ring r) {
        intvec * ww = to_intvec(w);
        ideal id = id_HomogenizeW(a, v, ww, r);
        return id;
    });

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

    Singular.method("id_MultSect", [](void * ids, int len, ring r) {
        const ring origin = currRing;
        rChangeCurrRing(r);
        ideal id = idMultSect(reinterpret_cast<resolvente>(ids), len);
        rChangeCurrRing(origin);
        return id;
    });

    Singular.method("id_Syzygies", &id_Syzygies_internal);

    Singular.method("id_sres", &id_sres_helper);

    Singular.method("id_fres", &id_fres_helper);

    Singular.method("id_res", &id_res_helper);

    Singular.method("id_Slimgb", &id_Slimgb_helper);

    Singular.method("id_MinStd", &id_MinStd_helper);
    Singular.method("id_TwoStd", &id_TwoStd_helper);
    Singular.method("id_Std", &id_Std_helper);
    Singular.method("id_StdHilb", &id_StdHilb_helper);
    Singular.method("id_StdHilbWeighted", &id_StdHilbWeighted_helper);

    Singular.method("id_InterRed", &id_InterRed_helper);

    Singular.method("id_Eliminate", [](ideal m, poly p, ring o) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal res = idElimination(m, p);
        rChangeCurrRing(origin);
        return res;
    });

    Singular.method("id_DivRem", [](ideal m, ideal sm, ring o) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal factors;
        ideal res = idDivRem(sm, m, factors, NULL);
        rChangeCurrRing(origin);
        return std::make_tuple(res, factors);
    });

    Singular.method("id_DivRem", [](ideal m, ideal sm, ring o, int flag) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal factors;
        ideal res = idDivRem(sm, m, factors, NULL,flag);
        rChangeCurrRing(origin);
        return std::make_tuple(res, factors);
    });

    Singular.method("id_DivRem_Unit", [](ideal m, ideal sm, ring o) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal factors;
        ideal unit;
        ideal res = idDivRem(sm, m, factors, &unit);
        rChangeCurrRing(origin);
        return std::make_tuple(res, factors, unit);
    });

    Singular.method("id_DivRem_Unit", [](ideal m, ideal sm, ring o, int flag) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal factors;
        ideal unit;
        ideal res = idDivRem(sm, m, factors, &unit, flag);
        rChangeCurrRing(origin);
        return std::make_tuple(res, factors, unit);
    });

    Singular.method("id_Lift", [](ideal m, ideal sm, ring o) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal rest;
        ideal res = idLift(m, sm, &rest, FALSE, FALSE);
        rChangeCurrRing(origin);
        return std::make_tuple(res, rest);
    });

    Singular.method("id_Lift", [](ideal m, ideal sm, bool goodShape,
                                  bool isSB, bool divide, ring o) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal rest;
        matrix unit;
        ideal res = idLift(m, sm, &rest, BOOLEAN(goodShape), BOOLEAN(isSB),
                                            BOOLEAN(divide), &unit, GbDefault);
        rChangeCurrRing(origin);
        return std::make_tuple(res, rest, unit);
    });

    Singular.method("id_LiftStd", [](ideal m, ring o, bool complete_reduction = false) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        matrix ma=mpNew(1,1);
        unsigned int crbit;
        if (complete_reduction)
            crbit = Sy_bit(OPT_REDSB);
        else
            crbit = 0;
        unsigned int save_opt = si_opt_1;
        si_opt_1 |= crbit;
        ideal res = idLiftStd(m, &ma, testHomog, NULL);
        si_opt_1 = save_opt;
        rChangeCurrRing(origin);
        return std::make_tuple(res, ma);
    });

    Singular.method("id_LiftStdSyz", [](ideal m, ring o, bool complete_reduction = false) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        matrix ma=mpNew(1,1);
        ideal syz=idInit(1,1);
        unsigned int crbit;
        if (complete_reduction)
            crbit = Sy_bit(OPT_REDSB);
        else
            crbit = 0;
        unsigned int save_opt = si_opt_1;
        si_opt_1 |= crbit;
        ideal res = idLiftStd(m, &ma, testHomog, &syz);
        si_opt_1 = save_opt;

        rChangeCurrRing(origin);
        return std::make_tuple(res, ma, syz);
    });

    Singular.method("id_Modulo", [](ideal a, ideal b, ring o) {
        const ring origin = currRing;
        rChangeCurrRing(o);
        ideal res = idModulo(a, b, testHomog);
        rChangeCurrRing(origin);
        return res;
    });

    Singular.method("id_Satstd", &id_Satstd);

    Singular.method("id_Saturation", [](ideal I, ideal J, ring r) {
        const ring origin = currRing;
        rChangeCurrRing(r);
        int d;
        ideal res = idSaturate(I, J, d, TRUE);
        rChangeCurrRing(origin);
        return std::make_tuple(res,d);
    });

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

    Singular.method("id_kbase", [](ideal I, int n, ring r) {
        ideal res;
        const ring origin = currRing;
        rChangeCurrRing(r);
        res = scKBase(n, I, r->qideal);
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
          ring im, void * cf_map) {
        const ring origin = currRing;
        rChangeCurrRing(pr);
        ideal I = maMapIdeal(map_id, pr, im_id, im, reinterpret_cast<nMapFunc>(cf_map));
        rChangeCurrRing(origin);
        return I;
    });
    Singular.method("idMinBase", [](ideal I, ring r) {
        const ring origin = currRing;
        rChangeCurrRing(r);
        ideal J = idMinBase(I);
        rChangeCurrRing(origin);
        return J;
    });
    Singular.method("scIndIndset", [](ideal I, ring r, jlcxx::ArrayRef<int> a, bool all) {
        const ring origin = currRing;
        rChangeCurrRing(r);
        lists L = scIndIndset(I, all, r->qideal);
        int n = rVar(r);
        int m = lSize(L);
        if(all == true && m >= 0)
        {
           for(int i = 0; i<=m; i++)
           {
              intvec * v = reinterpret_cast<intvec *>(L->m[i].data);
              int * content = v->ivGetVec();
              for(int j = 0; j < n; j++)
              {
                 a.push_back(content[j]);
              }
           }
        }
        else if(all == false && m >= 0)
        {
           intvec * v = reinterpret_cast<intvec *>(L->m[0].data);
           int * content = v->ivGetVec();
           for(int j = 0; j < n; j++)
           {
              a.push_back(content[j]);
           }
        }
        rChangeCurrRing(origin);
    });
    Singular.method("scDegree", [](ideal I, ring R)
    {
        const ring origin = currRing;
        rChangeCurrRing(R);
        SPrintStart();
        scDegree(I, NULL, R->qideal);
        char *s = SPrintEnd();
        s[strlen(s)-1]='\0';
        std::string res(s);
        omFree(s);
        rChangeCurrRing(origin);
        return res;
    });
    Singular.method("scDegree", [](ideal I, ring R, jlcxx::ArrayRef<int> w)
    {
        const ring origin = currRing;
        rChangeCurrRing(R);
        intvec * module_w = to_intvec(w);
        SPrintStart();
        scDegree(I, module_w, R->qideal);
        delete module_w;
        char *s = SPrintEnd();
        s[strlen(s)-1]='\0';
        std::string res(s);
        omFree(s);
        rChangeCurrRing(origin);
        return res;
    });
    Singular.method("scMultInt", [](ideal I, ring R) {
        const ring origin = currRing;
        rChangeCurrRing(R);
        int k = scMultInt(I, R->qideal);
        rChangeCurrRing(origin);
        return k;
    });
    Singular.method("scDimInt", [](ideal I, ring R) {
        const ring origin = currRing;
        rChangeCurrRing(R);
        int k = scDimInt(I, R->qideal);
        rChangeCurrRing(origin);
        return k;
    });
    Singular.method("scDimIntRing", [](ideal I, ring R) {
        const ring origin = currRing;
        rChangeCurrRing(R);
        int k = scDimIntRing(I, R->qideal);
        rChangeCurrRing(origin);
        return k;
    });
    Singular.method("fglmzero", [](ideal Isrc, ring Rsrc, ring Rdest ) {
        const ring origin = currRing;
        rChangeCurrRing(Rdest);
        ideal Idest = NULL;
        bool c = fglmzero(Rsrc, Isrc, Rdest, Idest, FALSE, FALSE);
        rChangeCurrRing(origin);
        return Idest;
    });
    Singular.method("scHilb", [](ideal I, ring r, jlcxx::ArrayRef<int> a) {
        const ring origin = currRing;
        rChangeCurrRing(r);
        intvec *v=hFirstSeries(I,NULL,r->qideal);
        int * content = v->ivGetVec();
        for(int j = 0; j < v->length(); j++)
          a.push_back(content[j]);
        delete v;
        rChangeCurrRing(origin);
    });
    Singular.method("scHilbWeighted", [](ideal I, ring r, jlcxx::ArrayRef<int> weights, jlcxx::ArrayRef<int> a) {
        intvec * w = to_intvec(weights);
        const ring origin = currRing;
        rChangeCurrRing(r);
        intvec *v=hFirstSeries(I,NULL,r->qideal,w);
        delete w;
        int * content = v->ivGetVec();
        for(int j = 0; j < v->length(); j++)
        {
          a.push_back(content[j]);
        }
        delete v;
        rChangeCurrRing(origin);
    });
    Singular.method("scHilbWeighted", [](ideal I, ring r, jlcxx::ArrayRef<int> weights, jlcxx::ArrayRef<int> shifts, jlcxx::ArrayRef<int> a) {
        intvec * w = to_intvec(weights);
        intvec * sh = to_intvec(shifts);
        const ring origin = currRing;
        rChangeCurrRing(r);
        intvec *v=hFirstSeries(I,sh,r->qideal,w);
        delete sh;
        delete w;
        int * content = v->ivGetVec();
        for(int j = 0; j < v->length(); j++)
        {
          a.push_back(content[j]);
        }
        delete v;
        rChangeCurrRing(origin);
    });
    Singular.method("scHilbPoly", [](ideal I, ring r, ring Qt) {
        const ring origin = currRing;
        rChangeCurrRing(r);
        poly h=hFirstSeries0p(I,r->qideal,w,r,Qt);
        rChangeCurrRing(origin);
	return h;
    });
    Singular.method("scHilbPolyWeighted", [](ideal I, ring r, jlcxx::ArrayRef<int> weights, ring Qt) {
        intvec * w = to_intvec(weights);
        const ring origin = currRing;
        rChangeCurrRing(r);
        poly h=hFirstSeries0p(I,r->qideal,w,r,Qt);
        delete w;
        rChangeCurrRing(origin);
	return h;
    });
    Singular.method("id_Homogen", id_Homogen);
    Singular.method("id_HomModule", [](jlcxx::ArrayRef<int> weights, ideal I, ring r) {
        intvec* w = NULL;
        bool res = id_HomModule(I, r->qideal, &w, r);
        if (w != NULL)
        {
            for (int i = 0; i < w->length(); i++)
                weights.push_back((*w)[i]);
            delete w;
        }
        return res;
    });
}
