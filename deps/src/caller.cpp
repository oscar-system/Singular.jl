#include "caller.h"

#include <Singular/tok.h>
#include <Singular/grammar.h>
#include <Singular/ipshell.h>
#include <Singular/lists.h>
#include <misc/intvec.h>

static jl_value_t * jl_int64_vector_type;
static jl_value_t * jl_int64_matrix_type;
static jl_value_t * jl_singular_number_type;
static jl_value_t * jl_singular_poly_type;
static jl_value_t * jl_singular_ring_type;
static jl_value_t * jl_singular_ideal_type;
static jl_value_t * jl_singular_matrix_type;
static jl_value_t * jl_singular_bigint_type;
static jl_value_t * jl_singular_bigintmat_type;
static jl_value_t * jl_singular_map_type;
static jl_value_t * jl_singular_resolution_type;
static jl_value_t * jl_singular_vector_type;

static jl_value_t * get_type_mapper()
{
    // clang-format off
    struct { int cmd; const char * name; } types[] = {
        {BIGINT_CMD, "BIGINT_CMD"},
        {NUMBER_CMD, "NUMBER_CMD"},
        {RING_CMD, "RING_CMD"},
        {POLY_CMD, "POLY_CMD"},
        {IDEAL_CMD, "IDEAL_CMD"},
        {MATRIX_CMD, "MATRIX_CMD"},
        {INT_CMD, "INT_CMD"},
        {STRING_CMD, "STRING_CMD"},
        {LIST_CMD, "LIST_CMD"},
        {INTMAT_CMD, "INTMAT_CMD"},
        {BIGINTMAT_CMD, "BIGINTMAT_CMD"},
        {MAP_CMD, "MAP_CMD"},
        {RESOLUTION_CMD, "RESOLUTION_CMD"},
        {MODUL_CMD, "MODUL_CMD"},
        {VECTOR_CMD, "VECTOR_CMD"},
        {INTVEC_CMD, "INTVEC_CMD"}};
    // clang-format on

    jl_array_t * return_array = jl_alloc_array_1d(
        jl_array_any_type, sizeof(types) / sizeof(types[0]));

    for (int i = 0; i < sizeof(types) / sizeof(types[0]); i++) {
        jl_array_t * current_return = jl_alloc_array_1d(jl_array_any_type, 2);
        jl_arrayset(current_return, jl_box_int64(types[i].cmd), 0);
        jl_arrayset(current_return,
                    reinterpret_cast<jl_value_t *>(jl_symbol(types[i].name)),
                    1);
        jl_arrayset(return_array,
                    reinterpret_cast<jl_value_t *>(current_return), i);
    }
    return reinterpret_cast<jl_value_t *>(return_array);
}

static void initialize_jl_c_types(jl_value_t * module_value)
{
    jl_module_t * module = reinterpret_cast<jl_module_t *>(module_value);
    jl_int64_vector_type =
        jl_apply_array_type((jl_value_t *)jl_int64_type, 1);
    jl_int64_matrix_type =
        jl_apply_array_type((jl_value_t *)jl_int64_type, 2);
    jl_singular_number_type = jl_get_global(module, jl_symbol("number"));
    jl_singular_poly_type = jl_get_global(module, jl_symbol("poly"));
    jl_singular_ring_type = jl_get_global(module, jl_symbol("ring"));
    jl_singular_ideal_type = jl_get_global(module, jl_symbol("ideal"));
    jl_singular_matrix_type = jl_get_global(module, jl_symbol("ip_smatrix"));
    jl_singular_bigint_type =
        jl_get_global(module, jl_symbol("__mpz_struct"));
    jl_singular_bigintmat_type =
        jl_get_global(module, jl_symbol("bigintmat"));
    jl_singular_map_type = jl_get_global(module, jl_symbol("sip_smap"));
    jl_singular_resolution_type =
        jl_get_global(module, jl_symbol("resolvente"));
}

static inline void * get_ptr_from_cxxwrap_obj(jl_value_t * obj)
{
    return *reinterpret_cast<void **>(obj);
}

// Safe way
// void* get_ptr_from_cxxwrap_obj(jl_value_t* obj){
//     return jl_unbox_voidpointer(jl_get_field(obj,"cpp_object"));
// }

jl_value_t * intvec_to_jl_array(intvec * v)
{
    int          size = v->length();
    jl_array_t * result = jl_alloc_array_1d(jl_int64_vector_type, size);
    int *        v_content = v->ivGetVec();
    for (int i = 0; i < size; i++) {
        jl_arrayset(result, jl_box_int64(static_cast<int64_t>(v_content[i])),
                    i);
    }
    return reinterpret_cast<jl_value_t *>(result);
}

jl_value_t * intmat_to_jl_array(intvec * v)
{
    int          rows = v->rows();
    int          cols = v->cols();
    jl_array_t * result = jl_alloc_array_2d(jl_int64_matrix_type, rows, cols);
    int64_t * result_ptr = reinterpret_cast<int64_t *> jl_array_data(result);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result_ptr[j + (i * cols)] = IMATELEM(*v, i, j);
        }
    }
    return reinterpret_cast<jl_value_t *>(result);
}

void * jl_array_to_intvec(jl_value_t * array_val)
{
    jl_array_t * array = reinterpret_cast<jl_array_t *>(array_val);
    int          size = jl_array_len(array);
    intvec *     result = new intvec(size);
    int *        result_content = result->ivGetVec();
    for (int i = 0; i < size; i++) {
        jl_value_t * current_entry = jl_arrayref(array, i);
        if (jl_is_int32(current_entry)) {
            result_content[i] =
                static_cast<int>(jl_unbox_int32(current_entry));
        }
        else if (jl_is_int64(current_entry)) {
            int64_t current_int64 = jl_unbox_int64(current_entry);
            result_content[i] = static_cast<int>(current_int64);
            if (result_content[i] != current_int64) {
                jl_error("Input entry does not fit in 32 bit integer");
            }
        }
    }
    return reinterpret_cast<void *>(result);
}

void * jl_array_to_intmat(jl_value_t * array_val)
{
    jl_array_t * array = reinterpret_cast<jl_array_t *>(array_val);
    int          rows = jl_array_dim(array, 0);
    int          cols = jl_array_dim(array, 1);
    intvec *     result = new intvec(rows, cols, 0);
    if (!jl_subtype(reinterpret_cast<jl_value_t *>(jl_typeof(array_val)),
                    reinterpret_cast<jl_value_t *>(jl_int64_matrix_type))) {
        jl_error("Input is not an Int64 matrix");
    }
    int64_t * array_data = reinterpret_cast<int64_t *>(jl_array_data(array));
    int *     vec_data = result->ivGetVec();
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            int64_t current_elem = array_data[j + (i * rows)];
            int     current_elem_32 = static_cast<int>(current_elem);
            if (current_elem != current_elem_32) {
                jl_error("Input entry does not fit in 32 bit integer");
            }
            IMATELEM(*result, i + 1, j + 1) = current_elem_32;
        }
    }
    return reinterpret_cast<void *>(result);
}

lists jl_array_to_list_helper(jl_value_t * args_val, jl_value_t * types_val)
{
    jl_array_t * args = reinterpret_cast<jl_array_t *>(args_val);
    jl_array_t * types = reinterpret_cast<jl_array_t *>(types_val);
    int          size = jl_array_len(args);
    lists      result = (lists)omAllocBin(slists_bin);
    result->Init(size);

    for (int i = 0; i < size; i++) {
       result->m[i].rtyp = static_cast<int>(jl_unbox_int64(jl_arrayref(types, i)));
       result->m[i].data = jl_unbox_voidpointer(jl_arrayref(args, i));
    }
    return result;
}

static void * get_ring_ref(ring r)
{
    // Since a call to a Singular library function destroys its arguments,
    // the call will decrease the number of references to the ring. So we
    // increase the reference count.
    r->ref++;
    return reinterpret_cast<void *>(r);
}

static jl_value_t * copy_polyptr_to_void(poly p, ring r)
{
    poly p_copy = p_Copy(p, r);
    return jl_box_voidpointer(reinterpret_cast<void *>(p_copy));
}

static jl_value_t * copy_idealptr_to_void(ideal i, ring r)
{
    ideal i_copy = id_Copy(i, r);
    return jl_box_voidpointer(reinterpret_cast<void *>(i_copy));
}

static jl_value_t * copy_bigintmatptr_to_void(bigintmat * m)
{
    bigintmat * m_copy = new bigintmat(m);
    return jl_box_voidpointer(reinterpret_cast<void *>(m_copy));
}

static void * copy_string_to_void(std::string s)
{
    return reinterpret_cast<void *>(omStrDup(s.c_str()));
}

bool translate_singular_type(jl_value_t * obj,
                             void **      args,
                             int *        argtypes,
                             int          i)
{
    jl_array_t * array = reinterpret_cast<jl_array_t *>(obj);
    int    cmd = static_cast<int>(jl_unbox_int64(jl_arrayref(array, 0)));
    void * arg = jl_unbox_voidpointer(jl_arrayref(array, 1));
    args[i] = arg;
    argtypes[i] = cmd;
    return true;
}

jl_value_t * get_julia_type_from_sleftv(leftv ret)
{
    jl_array_t * result = jl_alloc_array_1d(jl_array_any_type, 3);
    jl_arrayset(result, jl_false, 0);
    jl_arrayset(result, jl_box_voidpointer(ret->data), 1);
    ret->data = 0;
    jl_arrayset(result, jl_box_int64(ret->Typ()), 2);
    ret->rtyp = 0;
    return reinterpret_cast<jl_value_t *>(result);
}

jl_value_t * get_ring_content(ring r)
{
    ring save = currRing;
    rChangeCurrRing(r);

    // count elements
    idhdl h = r->idroot;
    int   nr = 0;
    while (h != NULL) {
        nr++;
        h = IDNEXT(h);
    }
    jl_array_t * result = jl_alloc_array_1d(jl_array_any_type, nr);
    h = r->idroot;
    nr = 0;
    while (h != NULL) {
        jl_array_t * current = jl_alloc_array_1d(jl_array_any_type, 3);
        jl_arrayset(current, jl_box_int64(IDTYP(h)), 0);
        jl_arrayset(current,
                    reinterpret_cast<jl_value_t *>(jl_symbol(IDID(h))), 1);
        {
            sleftv x; x.Copy((leftv)h);
            jl_arrayset(current, jl_box_voidpointer(x.data), 2);
        }
        jl_arrayset(result, reinterpret_cast<jl_value_t *>(current), nr);
        h = IDNEXT(h);
        nr++;
    }

    rChangeCurrRing(save);
    return reinterpret_cast<jl_value_t *>(result);
}

jl_value_t * call_singular_library_procedure(
    std::string s, ring r, jlcxx::ArrayRef<jl_value_t *> arguments)
{
    int    len = arguments.size();
    void * args[len];
    int    argtypes[len + 1];
    argtypes[len] = 0;
    for (int i = 0; i < len; i++) {
        bool result =
            translate_singular_type(arguments[i], args, argtypes, i);
        if (!result) {
            jl_error("Could not convert argument");
        }
    }
    BOOLEAN      err;
    jl_value_t * retObj;
    leftv        ret = ii_CallLibProcM(s.c_str(), args, argtypes, r, err);
    if (err) {
        jl_error("Could not call function");
    }
    if (ret->next != NULL) {
        int          len = ret->listLength();
        jl_array_t * list = jl_alloc_array_1d(jl_array_any_type, len + 1);
        jl_arrayset(list, jl_true, 0);
        for (int i = 0; i < len; ++i) {
            leftv next = ret->next;
            ret->next = 0;
            jl_arrayset(list, get_julia_type_from_sleftv(ret), i + 1);
            if (i > 0)
                omFreeBin(ret, sleftv_bin);
            ret = next;
        }
        retObj = reinterpret_cast<jl_value_t *>(list);
    }
    else {
        retObj = get_julia_type_from_sleftv(ret);
        omFreeBin(ret, sleftv_bin);
    }
    return retObj;
}

jl_value_t * call_singular_library_procedure_wo_rng(
    std::string name, void* rng, jlcxx::ArrayRef<jl_value_t *> arguments)
{
    return call_singular_library_procedure(name, reinterpret_cast<ring>(rng), arguments);
}

jl_value_t * lookup_singular_library_symbol_wo_rng(
    std::string pack,
    std::string name)
{
    int err = 2;
    jl_value_t * res = jl_nothing;
    jl_array_t * answer = jl_alloc_array_1d(jl_array_any_type, 2);
    leftv u = (leftv) IDROOT->get(pack.c_str(),0);
    if (u != NULL)
    {
        err--;
        idhdl v = ((package)(u->Data()))->idroot->get(name.c_str(), 0);
        if (v != NULL)
        {
            err--;
            sleftv x;
            x.Copy((leftv)v);
            res = get_julia_type_from_sleftv(&x);
        }
    }
    // return to julia [err, res]
    // err=0: no error, res is the value of the symbol
    // err=1: package found but symbol not found, res is junk
    // err=2: package not found, res is junk
    jl_arrayset(answer, jl_box_int64(err), 0);
    jl_arrayset(answer, res, 1);
    return reinterpret_cast<jl_value_t *>(answer);
}

jl_value_t * convert_nested_list(void * l_void)
{
    lists        l = reinterpret_cast<lists>(l_void);
    int          len = lSize(l) + 1;
    jl_array_t * result_array = jl_alloc_array_1d(jl_array_any_type, len);
    for (int i = 0; i < len; i++) {
        leftv current = &(l->m[i]);
        if (current->Typ() == LIST_CMD) {
            jl_arrayset(
                result_array,
                convert_nested_list(reinterpret_cast<void *>(current->data)),
                i);
        }
        else {
            jl_arrayset(result_array, get_julia_type_from_sleftv(current), i);
        }
    }
    return reinterpret_cast<jl_value_t *>(result_array);
}

void * create_syStrategy_data(syStrategy res, ring o)
{
    const ring origin = currRing;
    rChangeCurrRing(o);
    syStrategy temp = syCopy(res);
    rChangeCurrRing(origin);
    return reinterpret_cast<void *>(temp);
}

void singular_define_caller(jlcxx::Module & Singular)
{
    Singular.method("load_library", [](std::string name) {
        char * plib = iiConvName(name.c_str());
        idhdl  h = ggetid(plib);
        omFree(plib);
        if (h == NULL) {
            BOOLEAN bo = iiLibCmd(omStrDup(name.c_str()), TRUE, TRUE, FALSE);
            if (bo)
                return jl_false;
        }
        return jl_true;
    });
    Singular.method("lookup_singular_library_symbol_wo_rng",
                    &lookup_singular_library_symbol_wo_rng);
    Singular.method("call_singular_library_procedure",
                    &call_singular_library_procedure);
    Singular.method("call_singular_library_procedure",
                    &call_singular_library_procedure_wo_rng);
    Singular.method("get_type_mapper", &get_type_mapper);
    Singular.method("initialize_jl_c_types", &initialize_jl_c_types);


    Singular.method("NUMBER_CMD_CASTER",
                    [](void * obj) { return reinterpret_cast<number>(obj); });
    Singular.method("RING_CMD_CASTER",
                    [](void * obj) { return reinterpret_cast<ring>(obj); });
    Singular.method("POLY_CMD_CASTER",
                    [](void * obj) { return reinterpret_cast<poly>(obj); });
    Singular.method("IDEAL_CMD_CASTER",
                    [](void * obj) { return reinterpret_cast<ideal>(obj); });
    Singular.method("MATRIX_CMD_CASTER",
                    [](void * obj) { return reinterpret_cast<matrix>(obj); });
    Singular.method("INT_CMD_CASTER", [](void * obj) {
        return jl_box_int64(reinterpret_cast<long>(obj));
    });
    Singular.method("STRING_CMD_CASTER", [](void * obj) {
        return std::string(reinterpret_cast<char *>(obj));
    });
    Singular.method("INTVEC_CMD_CASTER", [](void * obj) {
        return intvec_to_jl_array(reinterpret_cast<intvec *>(obj));
    });
    Singular.method("INTMAT_CMD_CASTER", [](void * obj) {
        return intmat_to_jl_array(reinterpret_cast<intvec *>(obj));
    });
    Singular.method("BIGINT_CMD_CASTER", [](void * obj) {
        return reinterpret_cast<__mpz_struct *>(obj);
    });
    Singular.method("BIGINTMAT_CMD_CASTER", [](void * obj) {
        return reinterpret_cast<bigintmat *>(obj);
    });
    Singular.method("MAP_CMD_CASTER", [](void * obj) {
        return reinterpret_cast<sip_smap *>(obj);
    });
    Singular.method("RESOLUTION_CMD_CASTER", [](void * obj) {
        return reinterpret_cast<syStrategy>(obj);
    });
    Singular.method("LIST_CMD_TRAVERSAL", &convert_nested_list);
    Singular.method("get_ring_content", &get_ring_content);
    Singular.method("get_ring_ref", &get_ring_ref);
    Singular.method("copy_polyptr_to_void", &copy_polyptr_to_void);
    Singular.method("copy_idealptr_to_void", &copy_idealptr_to_void);
    Singular.method("copy_bigintmatptr_to_void", &copy_bigintmatptr_to_void);
    Singular.method("jl_array_to_intvec", &jl_array_to_intvec);
    Singular.method("jl_array_to_intmat", &jl_array_to_intmat);
    Singular.method("copy_string_to_void", &copy_string_to_void);
    Singular.method("jl_array_to_void", [] (jl_value_t * args_val,
                                                   jl_value_t * types_val,
                                                   ring R) {
        auto origin = currRing;
        rChangeCurrRing(R);
        lists l = jl_array_to_list_helper(args_val, types_val);
        rChangeCurrRing(origin);
        return (void *) l;
    });

    Singular.method("create_syStrategy_data", &create_syStrategy_data);

}
