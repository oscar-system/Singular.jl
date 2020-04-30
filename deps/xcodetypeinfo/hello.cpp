#include <string>
#include <vector>

#include "jlcxx/jlcxx.hpp"

std::string greet(const std::string& str)
{
   return "hello, world";
}

JLCXX_MODULE define_module_hello(jlcxx::Module& mod)
{
  mod.method("greet", &greet);
}
