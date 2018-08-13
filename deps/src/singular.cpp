#include <string>
#include <iostream>

#include "jlcxx/jlcxx.hpp"

#include <Singular/libsingular.h>
#include <Singular/cntrlc.h>


JULIA_CPP_MODULE_BEGIN(registry)
  jlcxx::Module& Singular = registry.create_module("libSingular");

  Singular.method("siInit",[](const char* path){ siInit(const_cast<char*>(path)); });
  Singular.method("versionString",[](){ return const_cast<const char*>(versionString()); });

  Singular.add_type<n_Procs_s>("coeffs");
//   Singular.add_type<snumber>("number");
JULIA_CPP_MODULE_END
