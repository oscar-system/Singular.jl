#include <string>
#include <iostream>

#include "jlcxx/jlcxx.hpp"

#include <Singular/libsingular.h>
#include <Singular/cntrlc.h>


void test_init(const char* test){
    siInit(const_cast<char*>(test));
}

const char* test_versionString(){
    return versionString();
}

JULIA_CPP_MODULE_BEGIN(registry)
  jlcxx::Module& Singular = registry.create_module("libSingular");
  Singular.method("si_init",&test_init);
  Singular.method("versionString",&test_versionString);

JULIA_CPP_MODULE_END
