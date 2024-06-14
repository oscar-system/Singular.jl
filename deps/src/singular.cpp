#include "includes.h"
#include "coeffs.h"
#include "rings.h"
#include "ideals.h"
#include "matrices.h"
#include "caller.h"
#include "coeff_rings.h"

static std::string              singular_return;
static std::string              singular_error;
static std::string              singular_warning;
std::vector<std::string> singular_error_log;

// Internal singular interpreter variable
extern int inerror;

// these are the temporary callbacks for calls to the interpreter
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

/*
   This is the non-temporary callback for all errors (unless the temporary
   ones are in use by call_interpreter). We would like to simultaneously:
    1. be able to check and report errors via libSingular.check_error()
    2. know when errors have been generated but uncaught by the julia code so
       that libSingular.check_error() can be inserted into the right place
   Unfortunately, a single call to the Singular kernel can generate multiple
   calls to WerrorS_callback, thus we don't know if previous errors were
   generated as a result of a missing libSingular.check_error() or if Singular
   has just called WerrorS_callback 10 times in the same function.
   The compromise here is to keep the full backlog of unreported errors and
   start complaining to stderr once the backlog gets too long.
*/
static void WerrorS_and_reset(const char * s)
{
  if (singular_error_log.size() > 9)
  {
    for (auto & si : singular_error_log)
      std::cerr << si << std::endl;
    std::cerr << "!!! Singular error(s) unhandled by julia !!!" << std::endl << std::endl;
  }
  singular_error_log.emplace_back(s);
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
  Singular.set_const("n_transExt", n_transExt);
  Singular.set_const("n_Nemo_AnticNumberField", n_Nemo_AnticNumberField);
  Singular.set_const("n_Nemo_QQField", n_Nemo_QQField);
  Singular.set_const("n_Nemo_ZZRing", n_Nemo_ZZRing);
  Singular.set_const("n_Nemo_FqPolyRepField", n_Nemo_FqPolyRepField);
  Singular.set_const("n_Nemo_fqPolyRepField", n_Nemo_fqPolyRepField);
  Singular.set_const("n_Nemo_Field", n_Nemo_Field);
  Singular.set_const("n_Nemo_Ring", n_Nemo_Ring);
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
  Singular.add_type<bigintmat>("bigintmat");

  /* monomial orderings */
  Singular.set_const("ringorder_no", ringorder_no);
  Singular.set_const("ringorder_lp", ringorder_lp);
  Singular.set_const("ringorder_ip", ringorder_ip);
  Singular.set_const("ringorder_dp", ringorder_dp);
  Singular.set_const("ringorder_Dp", ringorder_Dp);
  Singular.set_const("ringorder_wp", ringorder_wp);
  Singular.set_const("ringorder_Wp", ringorder_Wp);
  Singular.set_const("ringorder_Ip", ringorder_Ip);
  Singular.set_const("ringorder_ls", ringorder_ls);
  Singular.set_const("ringorder_is", ringorder_is);
  Singular.set_const("ringorder_ds", ringorder_ds);
  Singular.set_const("ringorder_Ds", ringorder_Ds);
  Singular.set_const("ringorder_ws", ringorder_ws);
  Singular.set_const("ringorder_Ws", ringorder_Ws);
  Singular.set_const("ringorder_a", ringorder_a);
  Singular.set_const("ringorder_M", ringorder_M);
  Singular.set_const("ringorder_c", ringorder_c);
  Singular.set_const("ringorder_C", ringorder_C);
  Singular.set_const("ringorder_s", ringorder_s);
  Singular.set_const("ringorder_S", ringorder_S);
  Singular.set_const("ringorder_IS", ringorder_IS);

  Singular.method("ringorder_to_int", [](rRingOrder_t a) {
    return static_cast<int>(a);
  });
  Singular.method("ringorder_from_int", [](int a) {
    return static_cast<rRingOrder_t>(a);
  });

  Singular.method("siInit", [](const char * path) {
    siInit(const_cast<char *>(path));
    WerrorS_callback = WerrorS_and_reset;
  });
  Singular.method("versionString", []() {
    return const_cast<const char *>(versionString());
  });
  Singular.method("version", []() {
    return SINGULAR_VERSION;
  });

  Singular.method("have_error", []() {
    return !singular_error_log.empty();
  });

  Singular.method("get_and_clear_error", []() {
    errorreported = 0;
    inerror = 0;
    std::stringstream ss;
    for (auto & si : singular_error_log)
      ss << si << std::endl;
    singular_error_log.clear();
    return ss.str();
  });

#define SETTER(A, B)                                                                     \
  else if (opt == #B)                                                                    \
  {                                                                                      \
    old_value = (A & Sy_bit(B)) != 0;                                                    \
    A = value ? (A | Sy_bit(B)) : (A & ~Sy_bit(B));                                      \
  }

  // all of the global setters return the previous value
  Singular.method("set_option", [](std::string opt, bool value) {
    bool old_value = false;
    if (false)
      ;
    SETTER(si_opt_2, V_QUIET)
    SETTER(si_opt_2, V_QRING)
    SETTER(si_opt_2, V_SHOW_MEM)
    SETTER(si_opt_2, V_YACC)
    SETTER(si_opt_2, V_REDEFINE)
    SETTER(si_opt_2, V_LOAD_LIB)
    SETTER(si_opt_2, V_DEBUG_LIB)
    SETTER(si_opt_2, V_LOAD_PROC)
    SETTER(si_opt_2, V_DEF_RES)
    SETTER(si_opt_2, V_SHOW_USE)
    SETTER(si_opt_2, V_IMAP)
    SETTER(si_opt_2, V_PROMPT)
    SETTER(si_opt_2, V_NSB)
    SETTER(si_opt_2, V_CONTENTSB)
    SETTER(si_opt_2, V_CANCELUNIT)
    SETTER(si_opt_2, V_MODPSOLVSB)
    SETTER(si_opt_2, V_UPTORADICAL)
    SETTER(si_opt_2, V_FINDMONOM)
    SETTER(si_opt_2, V_COEFSTRAT)
    SETTER(si_opt_2, V_IDLIFT)
    SETTER(si_opt_2, V_LENGTH)
    SETTER(si_opt_2, V_ALLWARN)
    SETTER(si_opt_2, V_INTERSECT_ELIM)
    SETTER(si_opt_2, V_INTERSECT_SYZ)
    SETTER(si_opt_2, V_DEG_STOP)

    SETTER(si_opt_1, OPT_PROT)
    SETTER(si_opt_1, OPT_REDSB)
    SETTER(si_opt_1, OPT_NOT_BUCKETS)
    SETTER(si_opt_1, OPT_NOT_SUGAR)
    SETTER(si_opt_1, OPT_INTERRUPT)
    SETTER(si_opt_1, OPT_SUGARCRIT)
    SETTER(si_opt_1, OPT_DEBUG)
    SETTER(si_opt_1, OPT_REDTHROUGH)
    SETTER(si_opt_1, OPT_NO_SYZ_MINIM)
    SETTER(si_opt_1, OPT_RETURN_SB)
    SETTER(si_opt_1, OPT_FASTHC)
    SETTER(si_opt_1, OPT_OLDSTD)
    SETTER(si_opt_1, OPT_STAIRCASEBOUND)
    SETTER(si_opt_1, OPT_MULTBOUND)
    SETTER(si_opt_1, OPT_DEGBOUND)
    SETTER(si_opt_1, OPT_REDTAIL)
    SETTER(si_opt_1, OPT_INTSTRATEGY)
    SETTER(si_opt_1, OPT_FINDET)
    SETTER(si_opt_1, OPT_INFREDTAIL)
    SETTER(si_opt_1, OPT_SB_1)
    SETTER(si_opt_1, OPT_NOTREGULARITY)
    SETTER(si_opt_1, OPT_WEIGHTM)
    else
    {
      std::cerr << "unknown option " << opt << std::endl;
    }
    return old_value;
  });
  // all of the global setters return the previous value
  Singular.method("set_option", [](std::string opt, bool value, ring r) {
    bool old_value = false;
    ring oldring=currRing;
    if (r!=NULL) rChangeCurrRing(r);
    if (false)
      ;
    SETTER(si_opt_2, V_QUIET)
    SETTER(si_opt_2, V_QRING)
    SETTER(si_opt_2, V_SHOW_MEM)
    SETTER(si_opt_2, V_YACC)
    SETTER(si_opt_2, V_REDEFINE)
    SETTER(si_opt_2, V_LOAD_LIB)
    SETTER(si_opt_2, V_DEBUG_LIB)
    SETTER(si_opt_2, V_LOAD_PROC)
    SETTER(si_opt_2, V_DEF_RES)
    SETTER(si_opt_2, V_SHOW_USE)
    SETTER(si_opt_2, V_IMAP)
    SETTER(si_opt_2, V_PROMPT)
    SETTER(si_opt_2, V_NSB)
    SETTER(si_opt_2, V_CONTENTSB)
    SETTER(si_opt_2, V_CANCELUNIT)
    SETTER(si_opt_2, V_MODPSOLVSB)
    SETTER(si_opt_2, V_UPTORADICAL)
    SETTER(si_opt_2, V_FINDMONOM)
    SETTER(si_opt_2, V_COEFSTRAT)
    SETTER(si_opt_2, V_IDLIFT)
    SETTER(si_opt_2, V_LENGTH)
    SETTER(si_opt_2, V_ALLWARN)
    SETTER(si_opt_2, V_INTERSECT_ELIM)
    SETTER(si_opt_2, V_INTERSECT_SYZ)
    SETTER(si_opt_2, V_DEG_STOP)

    SETTER(si_opt_1, OPT_PROT)
    SETTER(si_opt_1, OPT_REDSB)
    SETTER(si_opt_1, OPT_NOT_BUCKETS)
    SETTER(si_opt_1, OPT_NOT_SUGAR)
    SETTER(si_opt_1, OPT_INTERRUPT)
    SETTER(si_opt_1, OPT_SUGARCRIT)
    SETTER(si_opt_1, OPT_DEBUG)
    SETTER(si_opt_1, OPT_REDTHROUGH)
    SETTER(si_opt_1, OPT_NO_SYZ_MINIM)
    SETTER(si_opt_1, OPT_RETURN_SB)
    SETTER(si_opt_1, OPT_FASTHC)
    SETTER(si_opt_1, OPT_OLDSTD)
    SETTER(si_opt_1, OPT_STAIRCASEBOUND)
    SETTER(si_opt_1, OPT_MULTBOUND)
    SETTER(si_opt_1, OPT_DEGBOUND)
    SETTER(si_opt_1, OPT_REDTAIL)
    SETTER(si_opt_1, OPT_INTSTRATEGY)
    SETTER(si_opt_1, OPT_FINDET)
    SETTER(si_opt_1, OPT_INFREDTAIL)
    SETTER(si_opt_1, OPT_SB_1)
    SETTER(si_opt_1, OPT_NOTREGULARITY)
    SETTER(si_opt_1, OPT_WEIGHTM)
    else
    {
      std::cerr << "unknown option " << opt << std::endl;
    }
    if (r!=NULL)
    {
      r->options=si_opt_1;
      rChangeCurrRing(oldring);
    }
    return old_value;
  });

#undef SETTER

  // the "printlevel" system variable in Singular
  Singular.method("set_printlevel", [](int level) {
    int old_level = printlevel;
    printlevel = level;
    return old_level;
  });

  // the "degBound" system variable in Singular
  Singular.method("set_degBound", [](int degb) {
    int old_degb = Kstd1_deg;
    Kstd1_deg = degb;
    if (Kstd1_deg != 0)
      si_opt_1 |= Sy_bit(OPT_DEGBOUND);
    else
      si_opt_1 &= ~Sy_bit(OPT_DEGBOUND);
    return old_degb;
  });

  // the "multBound" system variable in Singular
  Singular.method("set_multBound", [](int mu) {
    int old_mu = Kstd1_mu;
    Kstd1_mu = mu;
    if (Kstd1_mu != 0)
      si_opt_1 |= Sy_bit(OPT_MULTBOUND);
    else
      si_opt_1 &= ~Sy_bit(OPT_MULTBOUND);
    return old_mu;
  });

  Singular.method("set_randomseed", [](int m) {
    int old_m = siSeed;
    if (m != 0)
    {
      siSeed = m;
      factoryseed(m);
    }
    return old_m;
  });

  Singular.method("random", []() {
    return siRand();
  });

  singular_define_coeffs(Singular);
  singular_define_rings(Singular);
  singular_define_ideals(Singular);
  singular_define_matrices(Singular);
  singular_define_caller(Singular);
  singular_define_coeff_rings(Singular);


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
    bool        err = iiAllStart(NULL, const_cast<char *>(input_str.c_str()), BT_proc, 0);
    inerror = 0;
    errorreported = 0;

    // get output
    jl_array_t * result = jl_alloc_array_1d(jl_array_any_type, 4);
    jl_array_ptr_set(result, 0, err ? jl_true : jl_false);
    jl_array_ptr_set(result, 1, jl_cstr_to_string(singular_return.c_str()));
    jl_array_ptr_set(result, 2, jl_cstr_to_string(singular_error.c_str()));
    jl_array_ptr_set(result, 3, jl_cstr_to_string(singular_warning.c_str()));

    // restore old callbacks
    PrintS_callback = default_print;
    WerrorS_callback = default_error;
    WarnS_callback = default_warning;

    return reinterpret_cast<jl_value_t *>(result);
  });

  /****************************
   ** from resolutions.jl
   ***************************/

  Singular.method("res_Delete_helper", [](syStrategy ra, ring o) {
    syKillComputation(ra, o);
  });

  Singular.method("res_Copy", [](syStrategy ra, ring o) {
    const ring origin = currRing;
    rChangeCurrRing_wo_options(o);
    syStrategy temp = syCopy(ra);
    rChangeCurrRing_wo_options(origin);
    return temp;
  });

  Singular.method("getindex_internal", [](syStrategy ra, int64_t k, bool minimal) {
    if (minimal)
    {
      return ra->minres[k];
    }
    return (ideal)ra->fullres[k];
  });

  Singular.method("syMinimize", [](syStrategy ra, ring o) {
    const ring origin = currRing;
    rChangeCurrRing_wo_options(o);
    syStrategy result = syCopy(ra);
    syMinimize(result);
    rChangeCurrRing_wo_options(origin);
    return result;
  });

  Singular.method("syMinimize_map", [](syStrategy ra, ring o) {
      const ring origin = currRing;
      rChangeCurrRing_wo_options(o);
      ideal T=NULL;
      syMinimize_with_map(ra,T);
      matrix TT=id_Module2Matrix(T,o);
      rChangeCurrRing_wo_options(origin);
      return TT;
  });

  Singular.method("get_minimal_res", [](syStrategy ra) {
        return reinterpret_cast<void *>(ra->minres);
  });

  Singular.method("get_full_res", [](syStrategy ra) {
    return reinterpret_cast<void *>(ra->fullres);
  });

  Singular.method("get_sySize", [](syStrategy ra) {
    return static_cast<int64_t>(sySize(ra));
  });

  Singular.method("create_SyStrategy", [](void * res_void, int64_t len, ring r) {
    resolvente res = reinterpret_cast<resolvente>(res_void);
    syStrategy result = (syStrategy)omAlloc0(sizeof(ssyStrategy));
    result->list_length = static_cast<short>(len);
    result->length = static_cast<int>(len);
    resolvente res_cp = (resolvente)omAlloc0((len + 1) * sizeof(ideal));
    for (int i = 0; i < len; i++)
    {
      if (res[i] != NULL)
      {
        res_cp[i] = id_Copy(res[i], r);
      }
    }
    result->fullres = res_cp;
    result->syRing = r;
    return result;
  });

  Singular.method("syBetti_internal", [](void * ra, int len, ring o) {
    const ring origin = currRing;
    rChangeCurrRing_wo_options(o);
    int      dummy;
    intvec * iv =
        syBetti(reinterpret_cast<resolvente>(ra), len, &dummy, NULL, FALSE, NULL);
    rChangeCurrRing_wo_options(origin);
    int  nrows = iv->rows();
    int  ncols = iv->cols();
    auto betti = (int *)malloc(ncols * nrows * sizeof(int));
    for (int i = 0; i < ncols; i++)
    {
      for (int j = 0; j < nrows; j++)
      {
        betti[i * nrows + j] = IMATELEM(*iv, j + 1, i + 1);
      }
    }
    delete (iv);
    return std::make_tuple(betti, nrows, ncols);
  });

  Singular.method("PrintS", &PrintS);
  Singular.method("StringAppendS", &StringAppendS);
}
