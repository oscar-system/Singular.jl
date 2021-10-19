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

/*
    calling rDefault with a non-null wvhdl argument:
    All ordering data should be passed in a serialized ord_data array. It is
    then deserialized into the rDefault arguments.
*/
auto rDefault_wvhdl_helper(coeffs cf,
                           jlcxx::ArrayRef<uint8_t *> vars,
                           jlcxx::ArrayRef<int> ord_data,
                           unsigned long bitmask)
{
    int nvars = vars.size();
    char ** vars_ptr = new char *[nvars];

    for (int i = 0; i < nvars; i++)
        vars_ptr[i] = reinterpret_cast<char *>(vars[i]);

    size_t j = 0;
    int nord = ord_data[j++];

    // rDefault is going to take ownership of all 4 of these pointers
    rRingOrder_t * ord = (rRingOrder_t *) omAlloc0((nord + 1) * sizeof(rRingOrder_t));
    int * blk0 = (int *) omAlloc0((nord + 1) * sizeof(int));
    int * blk1 = (int *) omAlloc0((nord + 1) * sizeof(int));
    int ** wvhdl = (int **) omAlloc0((nord + 1) * sizeof(int *));

    for (int i = 0; i < nord; i++)
    {
        ord[i] = (rRingOrder_t) ord_data[j++];
        blk0[i] = ord_data[j++];
        blk1[i] = ord_data[j++];
        int len = ord_data[j++];
        if (len > 0)
        {
            wvhdl[i] = (int *) omAlloc0(len * sizeof(int));
            for (int k = 0; k < len; k++)
                wvhdl[i][k] = ord_data[j++];
        }
        else
        {
            wvhdl[i] = nullptr;
        }
    }

    auto r = rDefault(cf, nvars, vars_ptr, nord, ord, blk0, blk1, wvhdl, bitmask);
    delete[] vars_ptr;
    r->ShortOut = 0;

    return r;
}

/*
    Serialize the ord, blk0, blk1, wvhdl data
    call with ord_data.size() == 0 since ArrayRef doesn't support resize
*/
void rOrdering_helper(jlcxx::ArrayRef<int> ord_data, const ring r)
{
    int l;
    size_t initial_size = ord_data.size();

    ord_data.push_back(0);

    for (l = 0; r->order[l]; l++)
    {
        rRingOrder_t ord = r->order[l];
        ord_data.push_back(static_cast<int>(ord));
        ord_data.push_back(r->block0[l]);
        ord_data.push_back(r->block1[l]);
        if (r->wvhdl[l] == nullptr)
        {
            ord_data.push_back(0);
        }
        else if (ord == ringorder_wp || ord == ringorder_Wp ||
                ord == ringorder_ws || ord == ringorder_Ws ||
                ord == ringorder_a || ord == ringorder_am || ord == ringorder_M)
        {
            int n = r->block1[l] - r->block0[l] + 1;
            if (ord == ringorder_M)
                n = n*n;
            ord_data.push_back(n);
            for (int i = 0; i < n; i++)
                ord_data.push_back(r->wvhdl[l][i]);
        }
        else
        {
            ord_data.push_back(0);
        }
    }

    ord_data[initial_size] = l;
}

auto rDefault_long_helper(coeffs                        cf,
                          jlcxx::ArrayRef<uint8_t *>    vars,
                          jlcxx::ArrayRef<rRingOrder_t> ord,
                          int *                         blk0,
                          int *                         blk1,
                          unsigned long                 bitmask)
{
    auto len = vars.size();
    char ** vars_ptr = new char *[len];
    for (int i = 0; i < len; i++) {
        vars_ptr[i] = reinterpret_cast<char *>(vars[i]);
        // std::strcpy(vars_ptr[i],vars[i].c_str());
    }
    auto len_ord = ord.size();
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

auto rDefault_Weyl_helper(coeffs                   cf,
                     jlcxx::ArrayRef<uint8_t *>    vars,
                     jlcxx::ArrayRef<rRingOrder_t> ord,
                     int *                         blk0,
                     int *                         blk1,
                     unsigned long                 bitmask)
{
    auto r = rDefault_long_helper(cf, vars, ord, blk0, blk1, bitmask);

    // vars = indeterminates x1, ..., xn and then partials Dx1, ..., Dxn
    int n = r->N/2;
    assume(2*n == r->N);
    matrix C = mpNew(2*n, 2*n);
    matrix D = mpNew(2*n, 2*n);

    for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++)
    {
        if (i < j)
        {
            // xj*xi = xi*xj
            MATELEM(C,i,j) = p_One(r);

            // Dxj*Dxi = Dxi*Dxj
            MATELEM(C,n+i,n+j) = p_One(r);
        }

        // Dxj*xi = xi*Dxj + bool(i == j)
        MATELEM(C,i,n+j) = p_One(r);
        if (i == j)
            MATELEM(D,i,n+j) = p_One(r);
    }

    ring R = rCopy(r);
    nc_CallPlural(C,D,NULL,NULL,R,true,false,true,r);

    // nc_CallPlural seems to take ownership of the matrices
    //mp_Delete(&C, r);
    //mp_Delete(&D, r);

    rDelete(r);
    return R;
}

// adapted from jiA_RING in ipassign.cc
static ring make_qring(ring r, ideal id)
{
    const ring origin = currRing;
    rChangeCurrRing(r);

    coeffs newcf = currRing->cf;

    const int cpos = idPosConstant(id);
    if (rField_is_Ring(currRing))
    {
        if (cpos >= 0)
        {
            newcf = n_CoeffRingQuot1(p_GetCoeff(id->m[cpos], currRing), currRing->cf);
            if (newcf == NULL)
            {
                rChangeCurrRing(origin);
                return NULL;
            }
        }
    }

    ring qr = rCopy(currRing);
    assume(qr->cf == currRing->cf);

    if (qr->cf != newcf)
    {
        nKillChar(qr->cf); // ???
        qr->cf = newcf;
    }

    ideal qid;

    if((rField_is_Ring(currRing)) && (cpos != -1))
    {
        int i, j;
        int *perm = (int *)omAlloc0((qr->N+1)*sizeof(int));

        for(i=qr->N;i>0;i--)
            perm[i]=i;

        nMapFunc nMap = n_SetMap(currRing->cf, newcf);
        qid = idInit(IDELEMS(id)-1,1);
        for (i = 0, j = 0; i<IDELEMS(id); i++)
            if (i != cpos)
                qid->m[j++] = p_PermPoly(id->m[i], perm, currRing, qr, nMap, NULL, 0);
    }
    else
    {
        qid = idrCopyR(id,currRing,qr);
    }

    idSkipZeroes(qid);

    if (currRing->qideal!=NULL) /* we are already in a qring! */
    {
        ideal tmp=idSimpleAdd(qid,currRing->qideal);
        // both ideals should be GB, so dSimpleAdd is sufficient
        idDelete(&qid);
        qid=tmp;
        // delete the qr copy of quotient ideal!!!
        idDelete(&qr->qideal);
    }
    if (idElem(qid)==0)
    {
        qr->qideal = NULL;
        id_Delete(&qid,currRing);
    }
    else
    {
        qr->qideal = qid;
    }

    // qr is a copy of currRing with the new qideal!
    #ifdef HAVE_PLURAL
    if (rIsPluralRing(currRing) &&(qr->qideal!=NULL))
    {
        if (nc_SetupQuotient(qr, currRing))
        {
            //WarnS("error in nc_SetupQuotient");
        }
    }
    #endif

    rChangeCurrRing(origin);
    return qr;
}


auto rDefault_Exterior_helper(coeffs               cf,
                     jlcxx::ArrayRef<uint8_t *>    vars,
                     jlcxx::ArrayRef<rRingOrder_t> ord,
                     int *                         blk0,
                     int *                         blk1,
                     unsigned long                 bitmask)
{
    auto r = rDefault_long_helper(cf, vars, ord, blk0, blk1, bitmask);
    int n = rVar(r);
    // first construct relations xj*xi = -xi*xj
    poly p = p_Neg(p_One(r), r);
    ring R = rCopy(r);
    nc_CallPlural(NULL,NULL,p,NULL,R,true,false,true,r);
    // then add xi^2 = 0
    ideal q = idInit(n,1);
    for (int i = 1; i <= n; i++)
        q->m[i-1] = p_Power(rGetVar(i, R), 2, R);
    ring S = make_qring(R, q);
    id_Delete(&q, R);
    p_Delete(&p, r);
    rDelete(R);
    rDelete(r);
    return S;
}

void singular_define_rings(jlcxx::Module & Singular)
{
    Singular.method("toPolyRef", [](void * ptr) {
       return reinterpret_cast<spolyrec*>(ptr);
    });
    Singular.method("freeAlgebra", &freeAlgebra);
    Singular.method("rDefault_helper", &rDefault_helper);
    Singular.method("rDefault_wvhdl_helper", &rDefault_wvhdl_helper);
    Singular.method("rOrdering_helper", &rOrdering_helper); // inverse of rDefault_wvhdl_helper
    Singular.method("rDefault_long_helper", &rDefault_long_helper);
    Singular.method("rDefault_Weyl_helper", &rDefault_Weyl_helper);
    Singular.method("rDefault_Exterior_helper", &rDefault_Exterior_helper);
    Singular.method("nc_CallPlural", [](matrix C, matrix D, ring r) {
        ring R = rCopy(r);
        nc_CallPlural(mp_Copy(C,r), mp_Copy(D,r), NULL, NULL, R, true, false, true, r);
        return R;
    });
    Singular.method("rDelete", &rDelete);
    Singular.method("rString", [](ip_sring * r) {
        auto s = rString(r);
        std::string ret_string(s);
        omFree(s);
        return ret_string;
    });
    Singular.method("rChar", &rChar);
    Singular.method("rGetVar", &rGetVar);
    Singular.method("rVar", &rVar);
    Singular.method("rIsLPRing", [](const ring r) {
        // return is 0 for not LP Ring, otherwise the number of variables for this ring
        #ifdef HAVE_SHIFTBBA
            return r->isLPring;
        #else
            return short(0);
        #endif
    });
    Singular.method("rRingVar", [](short i, const ring r) {
        return std::string(rRingVar(i, r));
    });
    Singular.method("rGetExpSize", [](unsigned long bitmask, int N) {
        int bits;
        return static_cast<unsigned int>(rGetExpSize(bitmask, bits, N));
    });
    Singular.method("rCoeffPtr", [](ring r){return r->cf;});
    Singular.method("rHasGlobalOrdering", &rHasGlobalOrdering);
    Singular.method("rHasMixedOrdering", &rHasMixedOrdering);
    Singular.method("rRing_ord_pure_dp", &rRing_ord_pure_dp);
    Singular.method("rRing_ord_pure_Dp", &rRing_ord_pure_Dp);
    Singular.method("rRing_ord_pure_lp", &rRing_ord_pure_lp);
    Singular.method("rIsQuotientRing", [](ring r) {
        return r->qideal != NULL;
    });
    Singular.method("rCopy", rCopy);
    Singular.method("rQuotientRing", [](ideal i, ring r) {
        ring Q = rCopy(r);
        Q->qideal = id_Copy(i, r);
        return Q;
    });
    Singular.method("rBitmask",
                    [](ip_sring * r) { return (unsigned int)r->bitmask; });
    Singular.method("rPar", [](coeffs cf){
                    return n_NumberOfParameters(cf);
    });
    Singular.method("p_Delete", [](spolyrec * p, ip_sring * r) {
        return p_Delete(&p, r);
    });
    Singular.method("p_Copy",
                    [](spolyrec * p, ip_sring * r) { return p_Copy(p, r); });
    Singular.method("p_IsOne", p_IsOne);
    Singular.method("p_One", p_One);
    Singular.method("p_IsUnit", p_IsUnit);
    Singular.method("p_GetExp", [](spolyrec * p, int i, ip_sring * r) {
        return p_GetExp(p, i, r);
    });
    Singular.method("p_GetComp", [](spolyrec * p, ip_sring * r) {
        return p_GetComp(p, r);
    });
    Singular.method("p_String", [](spolyrec * p, ip_sring * r) {
        auto s_ptr = p_String(p, r);
        std::string s(s_ptr);
        omFree(s_ptr);
        return s;
    });
    Singular.method("p_ISet", p_ISet);
    Singular.method("p_NSet", p_NSet);
    Singular.method("p_NSet",
                    [](void * p, ip_sring * r) { return p_NSet(reinterpret_cast<snumber*>(p), r); }
    );
    Singular.method("pLength", pLength);
    Singular.method("SetpNext",
                    [](spolyrec * p, spolyrec * q) { p->next = q; });
    Singular.method("pNext", [](spolyrec * a) {
        poly p = pNext(a);
        return p;
    });
    Singular.method("p_Init", [](ip_sring * r) { return p_Init(r); });
    Singular.method("p_Head", [](spolyrec * a, ip_sring * r) {
        poly p = p_Head(a, r); return p; });
    Singular.method("p_SetCoeff0", [](spolyrec * a, snumber * n, ip_sring * r) {
        p_SetCoeff0(a, n, r); });
    Singular.method("p_SetExp", [](spolyrec * a, int i, int v, ip_sring * r) {
        p_SetExp(a, i, v, r); });
    Singular.method("p_SetNext", [](spolyrec * a, spolyrec * m) {
        pNext(a) = m; });
    Singular.method("p_SortMerge", [](spolyrec * a, ip_sring * r) {
        return p_SortMerge(a, r); });
    Singular.method("p_SortAdd", p_SortAdd);
    Singular.method("p_Setm", p_Setm);
    Singular.method("p_Neg", p_Neg);
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
    Singular.method("p_Sub", p_Sub);
    Singular.method("p_Mult_q", p_Mult_q);
    Singular.method("pp_Mult_qq", pp_Mult_qq);
    Singular.method("p_Power", p_Power);
    Singular.method("p_Compare", p_Compare);
    Singular.method("p_EqualPolys",
                    [](spolyrec * p, spolyrec * q, ip_sring * r) {
                        return p_EqualPolys(p, q, r);
                    });
    Singular.method("p_Divide", p_Divide);
    Singular.method("p_DivRem", [](spolyrec * a, spolyrec * b, ip_sring * r) {
       poly rest;
       poly q = p_DivRem(a, b, rest, r);
       return std::make_tuple(reinterpret_cast<void *>(q), reinterpret_cast<void *>(rest));
    });
    Singular.method("p_Div_nn", p_Div_nn);
    Singular.method("p_Mult_nn", [](poly p, number n, ring r) {
       return p_Mult_nn(p, n, r);
    });
    Singular.method("p_LmIsConstant", [](poly p, ring r) {
       return p_LmIsConstant(p, r) != 0;
    });
    Singular.method("p_LmIsConstantComp", [](poly p, ring r) {
       return p_LmIsConstantComp(p, r) != 0;
    });
    Singular.method("p_IsDivisibleBy", [](spolyrec * p, spolyrec * q, ip_sring * r) {
       poly res;
       ideal I = idInit(1, 1);
       const ring origin = currRing;
       I->m[0] = q;
       rChangeCurrRing(r);
       res = kNF(I, NULL, p, 0, KSTD_NF_LAZY);
       rChangeCurrRing(origin);
       I->m[0] = NULL;
       id_Delete(&I, r);
       if (res == NULL)
          return true;
       else
       {
          p_Delete(&res, r);
          return false;
       }
    });
    Singular.method("singclap_gcd", singclap_gcd);

    Singular.method("singclap_extgcd", [](spolyrec * a, spolyrec * b, ip_sring * r) {
        poly res, s, t;
        singclap_extgcd(a, b, res, s, t, r);
        return std::make_tuple(res, s, t);
    });

    /* TODO check if we can get rid of this nonsensical function */
    Singular.method("p_ExtGcd_internal", [](spolyrec * a, spolyrec * b,
                                            void * res, void * s, void * t,
                                            ip_sring * r) {
        return singclap_extgcd(a, b, reinterpret_cast<spolyrec *&>(res),
                               reinterpret_cast<spolyrec *&>(s),
                               reinterpret_cast<spolyrec *&>(t), r);
    });
    Singular.method("singclap_sqrfree",
                    [](spolyrec * p, jlcxx::ArrayRef<int> a, ip_sring * r) {
                        const ring origin = currRing;
		        rChangeCurrRing(r);
			intvec * v = NULL;
			ideal I = singclap_sqrfree(pCopy(p), &v, 0, currRing);
			int * content = v->ivGetVec();
			for(int i=0; i<v->length(); i++)
			{
			  a.push_back(content[i]);
			}
			rChangeCurrRing(origin);
			delete v;
		        return I;
    });
    Singular.method("singclap_factorize",
                    [](spolyrec * p, jlcxx::ArrayRef<int> a, ip_sring * r) {
		    	const ring origin = currRing;
			rChangeCurrRing(r);
			intvec * v = NULL;
			ideal I = singclap_factorize(p_Copy(p,r), &v, 0, r);
			int * content = v->ivGetVec();
			for(int i=0; i<v->length(); i++)
			{
			  a.push_back(content[i]);
			}
			rChangeCurrRing(origin);
			delete v;
		        return I;
    });
    Singular.method("p_Content", p_Content);
    Singular.method("p_GetExpVL_internal", p_GetExpVL);
    Singular.method("p_GetExpVLV_internal", p_GetExpVLV);
    Singular.method("p_SetExpV_internal", p_SetExpV);
    Singular.method("p_SetExpVL_internal", p_SetExpVL);
    Singular.method("p_SetExpVLV_internal", p_SetExpVLV);
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

    Singular.method("p_Subst", [](poly p, int i, poly q, ring r) {
        poly p_cp = p_Copy(p, r);
        return p_Subst(p_cp, i, q, r);
    });
    Singular.method("maEvalAt", [](poly p, jlcxx::ArrayRef<snumber *> vals, ring r) {
       number * varr = (number *) omAlloc0(vals.size() * sizeof(number));
       for (int i = 0; i < vals.size(); i++)
          varr[i] = (number) vals[i];
       number res = maEvalAt(p, varr, r);
       omFree(varr);
       return res;
    });
    Singular.method("p_PermPoly", [](poly p, int * perm, ring old_ring,
                                     ring new_ring, void * map_func_ptr, int * par_perm) {
        nMapFunc map_func = reinterpret_cast<nMapFunc>(map_func_ptr);
        return p_PermPoly(p, perm, old_ring, new_ring, map_func, par_perm);
    });
    Singular.method("maFindPerm", [](ring src, jlcxx::ArrayRef<int> perm, ring dst,
                                     jlcxx::ArrayRef<int> par_perm){
        int *perm1 = (int *)omAlloc0((rVar(src)+1)*sizeof(int));
        int *par_perm1 = NULL;
        if (rPar(src) != 0) par_perm1 = (int *)omAlloc0((rPar(src) + 1)*sizeof(int));
        maFindPerm(src->names, rVar(src), rParameter(src), rPar(src),
                 dst->names, rVar(dst), rParameter(dst), rPar(dst),
                 perm1, par_perm1, dst->cf->type);
        for(int i = 0; i < rVar(src); i++)
        {
           perm.push_back(perm1[i]);
        }
        for(int j = 0; j < rPar(src); j++)
        {
           par_perm.push_back(par_perm1[j]);
        }
    });
    Singular.method("p_Jet",
                   [](poly p, int i, ring r) {
                       poly p_cp = p_Copy(p, r);
                       return p_Jet(p_cp, i, r);
    });
    Singular.method("p_Diff",
                   [](poly p, int i, ring r) {
                       poly p_cp = p_Copy(p, r);
                       return p_Diff(p_cp, i, r);
    });
    Singular.method("maMapPoly",
           [](poly map_p, ring pr, ideal im_id, ring im, void * cf_map) {
        return maMapPoly(map_p, pr, im_id, im, reinterpret_cast<nMapFunc>(cf_map));
    });
    Singular.method("p_GetOrder",
                   [](poly p, ring r) {
		       long res;
                       if( p != NULL)
		       {  res = p_GetOrder(p, r);}
		       else 
		       {  res = -1;}
		       return res; 
    });
}
