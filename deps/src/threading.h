#include <mutex>

#ifdef THREADSAFE_SINGULAR

namespace singularjl {
    extern std::recursive_mutex global_singular_lock;
}

#define ENTER_SINGULAR (singularjl::global_singular_lock.lock())
#define LEAVE_SINGULAR (singularjl::global_singular_lock.unlock())

#else

#define ENTER_SINGULAR ((void) 0)
#define LEAVE_SINGULAR ((void) 0)

#endif
