#include "matrices.h"

void singular_define_matrices(jlcxx::Module & Singular)
{
    Singular.method("ncols", [](matrix I) { return (int)MATCOLS(I); });

    Singular.method("nrows", [](matrix I) { return (int)MATROWS(I); });

    Singular.method("id_Module2Matrix", &id_Module2Matrix);

    Singular.method("id_Matrix2Module", &id_Matrix2Module);

    Singular.set_override_module(jl_base_module);
    Singular.method("getindex", [](matrix M, int i, int j) {
        return (poly)MATELEM(M, i, j);
    });
    Singular.unset_override_module();

    Singular.method("setindex", [](matrix M, poly p, int i, int j, ring R) {
        MATELEM(M, i, j) = p_Copy(p, R);
    });

    Singular.method("mp_Copy",
                    [](matrix M, ring R) { return mp_Copy(M, R); });

    Singular.method("mp_Delete",
                    [](matrix M, ring R) { return mp_Delete(&M, R); });

    Singular.method("mp_Add", &mp_Add);

    Singular.method("mp_Sub", &mp_Sub);

    Singular.method("mp_Transp", &mp_Transp);

    Singular.method("mp_Mult", &mp_Mult);

    Singular.method("mp_MultP", &mp_MultP);

    Singular.method("pMultMp", &pMultMp);

    Singular.method("mp_Equal", &mp_Equal);

    Singular.method("mpNew", [](int r, int c) {
        return mpNew(r, c);
    });

    Singular.method("mp_InitP", [](int n, poly p, ring R) {
        return mp_InitP(n, n, p_Copy(p, R), R);
    });

    Singular.method("mp_Wedge", [](matrix M,int n, ring R) {
        return mp_Wedge(M, n, R);
    });

    Singular.method("irrCharSeries", &singclap_irrCharSeries);

    Singular.method("iiStringMatrix", [](matrix I, int d, ring o) {
        auto str_ptr = iiStringMatrix(I, d, o);
        std::string s(iiStringMatrix(I, d, o));
        omFree(str_ptr);
        return s;
    });

    Singular.method("bigintmat_init", [](int r, int c) {
        return new bigintmat(r, c, coeffs_BIGINT);
    });
    Singular.method("bigintmat_clear", [](bigintmat * m) {
        delete m;
    });
    Singular.method("bigintmat_nrows", [](bigintmat * m) {
        return m->rows();
    });
    Singular.method("bigintmat_ncols", [](bigintmat * m) {
        return m->cols();
    });
    Singular.method("bigintmat_viewindex", [](bigintmat * m, int i, int j) {
        return m->view(i, j);
    });
    Singular.method("bigintmat_rawset", [](bigintmat * m, number n, int i, int j) {
        m->rawset(i, j, n, NULL);
    });
}
