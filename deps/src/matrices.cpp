#include "matrices.h"

void singular_define_matrices(jlcxx::Module & Singular)
{
    Singular.method("ncols", [](matrix I) { return (int)MATCOLS(I); });

    Singular.method("nrows", [](matrix I) { return (int)MATROWS(I); });

    Singular.method("id_Module2Matrix", &id_Module2Matrix);

    Singular.method("id_Matrix2Module", &id_Matrix2Module);

    Singular.method("getindex", [](matrix M, int i, int j) {
        return (poly)MATELEM(M, i, j);
    });

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

    Singular.method("mp_Equal", &mp_Equal);

    Singular.method("mpNew", [](ring R, int r, int c) {
        rChangeCurrRing(R);
        return mpNew(r, c);
    });

    Singular.method("mp_InitP", [](int n, poly p, ring R) {
        rChangeCurrRing(R);
        return mp_InitP(n, n, pCopy(p), R);
    });

    Singular.method("iiStringMatrix", [](matrix I, int d, ring o) {
        auto str_ptr = iiStringMatrix(I, d, o);
        std::string s(iiStringMatrix(I, d, o));
        omFree(str_ptr);
        return s;
    });
}
