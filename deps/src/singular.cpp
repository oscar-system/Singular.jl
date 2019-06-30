#include "jlcxx/jlcxx.hpp"
#include "includes.h"
#include "coeffs.h"
#include "rings.h"
#include "ideals.h"
#include "matrices.h"

static std::string singular_return;
static std::string singular_error;
static std::string singular_warning;

// Internal singular interpreter variable
extern int         inerror;

static void WerrorS_for_julia(const char * s)
{
    singular_error += s;
}

static void PrintS_for_julia(const char * s)
{
    singular_return += s;
}

static void WarningS_for_julia(const char * s)
{
    singular_warning += s;
}

JLCXX_MODULE define_julia_module(jlcxx::Module & Singular)
{
    Singular.add_type<n_Procs_s>("coeffs");
    Singular.add_bits<n_coeffType>("n_coeffType");
    Singular.set_const("n_Z", n_Z);
    Singular.set_const("n_Q", n_Q);
    Singular.set_const("n_Zn", n_Zn);
    Singular.set_const("n_Zp", n_Zp);
    Singular.set_const("n_GF", n_GF);
    Singular.set_const("n_unknown", n_unknown);
    Singular.add_type<snumber>("number");
    Singular.add_type<__mpz_struct>("__mpz_struct");
    Singular.add_type<ip_sring>("ring");
    Singular.add_type<spolyrec>("poly");
    // Singular.add_type<nMapFunc>("nMapFunc");
    // Singular.add_type<spolyrec>("vector");
    Singular.add_bits<rRingOrder_t>("rRingOrder_t");
    Singular.add_type<sip_sideal>("ideal");
    Singular.add_type<ip_smatrix>("ip_smatrix");
    Singular.add_type<ssyStrategy>("syStrategy");
    Singular.add_type<sip_smap>("sip_smap");

    /* monomial orderings */
    Singular.set_const("ringorder_no", ringorder_no);
    Singular.set_const("ringorder_lp", ringorder_lp);
    Singular.set_const("ringorder_rp", ringorder_rp);
    Singular.set_const("ringorder_dp", ringorder_dp);
    Singular.set_const("ringorder_Dp", ringorder_Dp);
    Singular.set_const("ringorder_ls", ringorder_ls);
    Singular.set_const("ringorder_rs", ringorder_rs);
    Singular.set_const("ringorder_ds", ringorder_ds);
    Singular.set_const("ringorder_Ds", ringorder_Ds);
    Singular.set_const("ringorder_c", ringorder_c);
    Singular.set_const("ringorder_C", ringorder_C);

    Singular.method("siInit", [](const char * path) {
        siInit(const_cast<char *>(path));
    });
    Singular.method("versionString", []() {
        return const_cast<const char *>(versionString());
    });

    singular_define_coeffs(Singular);
    singular_define_rings(Singular);
    singular_define_ideals(Singular);
    singular_define_matrices(Singular);


    // Calls the Singular interpreter with `input`.
    // `input` needs to be valid Singular input. 
    // Returns a 4-tuple:
    // 1. entry is a bool, indicated if an error has happened
    // 2. entry is the output as a string
    // 3. entry is the error output as a string
    // 4. entry is the warning output as a string
    Singular.method("call_interpreter", [](std::string input) {
        // save callbacks
        auto default_print = PrintS_callback;
        auto default_error = WerrorS_callback;
        auto default_warning = WarnS_callback;

        // set temporary new callbacks
        PrintS_callback = PrintS_for_julia;
        WerrorS_callback = WerrorS_for_julia;
        WarnS_callback = WarningS_for_julia;

        // cleanup return strings
        singular_return.clear();
        singular_error.clear();
        singular_warning.clear();

        // call interpreter
        std::string input_str = input + "\nreturn();";
        bool err = iiAllStart(NULL, const_cast<char *>(input_str.c_str()),
                              BT_proc, 0);
        inerror = 0;
        errorreported = 0;

        // get output
        jl_array_t * result = jl_alloc_array_1d(jl_array_any_type, 4);
        jl_arrayset(result, err ? jl_true : jl_false, 0);
        jl_arrayset(result, jl_cstr_to_string(singular_return.c_str()), 1);
        jl_arrayset(result, jl_cstr_to_string(singular_error.c_str()), 2);
        jl_arrayset(result, jl_cstr_to_string(singular_warning.c_str()), 3);

        // restore old callbacks
        PrintS_callback = default_print;
        WerrorS_callback = default_error;
        WarnS_callback = default_warning;

        return reinterpret_cast<jl_value_t *>(result);
    });

    /****************************
     ** from resolutions.jl
     ***************************/

    Singular.method("res_Delete_helper", [](void * ra_void, int len, ring o) {
        auto ra = reinterpret_cast<resolvente>(ra_void);
        for (int i = 0; i < len; i++) {
            id_Delete(ra + i, o);
        }
        omFreeSize((ADDRESS)ra, (len + 1) * sizeof(ideal));
    });

    Singular.method("res_Copy", [](void * ra_void, int len, ring o) {
        auto       ra = reinterpret_cast<resolvente>(ra_void);
        resolvente res = (resolvente)omAlloc0((len + 1) * sizeof(ideal));
        rChangeCurrRing(o);
        for (int i = len - 1; i >= 0; i--) {
            if (ra[i] != NULL)
                res[i] = id_Copy(ra[i], o);
        }
        return reinterpret_cast<void *>(res);
    });


    Singular.method("getindex", [](void * ra_void, int k) {
        auto ra = reinterpret_cast<resolvente>(ra_void);
        return (ideal)ra[k];
    });

    Singular.method("syMinimize", [](void * ra_void, int len, ring o) {
        auto       ra = reinterpret_cast<resolvente>(ra_void);
        const ring origin = currRing;
        syStrategy temp = (syStrategy)omAlloc0(sizeof(ssyStrategy));
        resolvente result;
        rChangeCurrRing(o);
        temp->fullres = (resolvente)omAlloc0((len + 1) * sizeof(ideal));
        for (int i = len - 1; i >= 0; i--) {
            if (ra[i] != NULL)
                temp->fullres[i] = idCopy(ra[i]);
        }
        temp->length = len;
        syMinimize(temp);
        result = temp->minres;
        temp->minres = NULL;
        // syMinimize increments this as it returns a value we ignore
        temp->references--;
        syKillComputation(temp, o);
        rChangeCurrRing(origin);
        return reinterpret_cast<void *>(result);
    });


    Singular.method("syBetti", [](void * ra_void, int len, ring o) {
        auto       ra = reinterpret_cast<resolvente>(ra_void);
        const ring origin = currRing;
        rChangeCurrRing(o);
        int      dummy;
        intvec * iv = syBetti(ra, len, &dummy, NULL, FALSE, NULL);
        rChangeCurrRing(origin);
        int  nrows = iv->rows();
        int  ncols = iv->cols();
        auto betti = (int *)malloc(ncols * nrows * sizeof(int));
        for (int i = 0; i < ncols; i++) {
            for (int j = 0; j < nrows; j++) {
                betti[i * nrows + j] = IMATELEM(*iv, j + 1, i + 1);
            }
        }
        delete (iv);
        return std::make_tuple(betti, nrows, ncols);
    });
}
