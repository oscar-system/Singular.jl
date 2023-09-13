#ifndef XCODE_JLCXX_WORKAROUND
#define XCODE_JLCXX_WORKAROUND

// This must be the very first include
#ifdef _LIBCPP_CONFIG
#error This header must be included before any system headers!
#endif

// Work around Xcode 11.4 issue until upstream libc++ fix arrives in Xcode:
// https://github.com/llvm/llvm-project/commit/2464d8135e
//
// First include __config from libc++, then override typeinfo flag
// to force use of address as hash instead of hashing the string.

#if defined(__APPLE__) && defined(FORCE_XCODE_TYPEINFO_MERGED)
#include <__config>
#if defined(_LIBCPP_HAS_MERGED_TYPEINFO_NAMES_DEFAULT) &&                                \
    _LIBCPP_HAS_MERGED_TYPEINFO_NAMES_DEFAULT == 0
#undef _LIBCPP_HAS_MERGED_TYPEINFO_NAMES_DEFAULT
#define _LIBCPP_HAS_MERGED_TYPEINFO_NAMES_DEFAULT 1
#else
#error Trying to work around Xcode 11.4 bug but libc++ macro not set as expected! \
 Please try rebuilding and create an issue if this reappears.
#endif
#endif

#include "jlcxx/jlcxx.hpp"

#endif
