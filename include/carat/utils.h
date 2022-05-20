#ifndef CARAT_UTILS_H
#define CARAT_UTILS_H

#include "config.h"

#include <stdlib.h>

#if defined(__GNUC__)
# define CARAT_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
#  define CARAT_CURRENT_FUNCTION	__func__
#elif defined(_MSC_VER)
#  define CARAT_CURRENT_FUNCTION __FUNCTION__
#else
#  define CARAT_CURRENT_FUNCTION "<unknown>"
#endif

extern void *xmalloc_(size_t size, const char *funcname);
extern void *xrealloc_(void *ptr, size_t size, const char *funcname);

#define xmalloc(size)  xmalloc_(size, CARAT_CURRENT_FUNCTION)
#define xrealloc(ptr, size)  xrealloc_(ptr, size, CARAT_CURRENT_FUNCTION)


#ifndef HAVE_STRLCPY
size_t strlcpy(char * dst, const char * src, size_t len);
#endif

#ifndef HAVE_STRLCAT
size_t strlcat(char * dst, const char * src, size_t len);
#endif


#endif /* CARAT_UTILS_H */
