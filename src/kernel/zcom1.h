/*
  common routines

  Copyright (c) 2006-2012 Cheng Zhang

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Usage:

  1.  It is designed for quick programming.
      For simple use, include this file and all functions will be available.
      But there might be many compiler warnings for unused functions.

  2.  You can include this file multiple times in a single file.

  3.  Function are static by default. To export functions,
      e.g., to make it easier to debug, or to avoid warnings of unused functions,
      define ZCOM_XFUNCS before the first inclusion.

  4.  To hand-pick specific set of modules, e.g.,
        #define ZCOM_PICK
        #define ZCOM_RNG
        #define ZCOM_ARGOPT
      before including this file, so other modules are skipped.

  5.  If the compiler supports keywords inline and restrict,
        #define INLINE inline
        #define RESRICT restrict
      before including this file.

  6.  Define HAVEVAM if the compiler supports variable-argument macros.

  7.  The def module defines `real' as a double, to override
        typedef float real;
        #define HAVEREAL 1
      before including this file (or equivalently define HAVE_REAL)
*/

/* ZCOM_PICK or ZCOM_NONE is used include only subset of modules to
 * 1. reduce the # of warnings for unused functions
 * 2. accelerate the compiling
 * 3. avoid multiple inclusions
 * By default, ZCOM_PICK is undefined, so everything is used. */
#ifdef ZCOM_NONE  /* equivalent to ZCOM_PICK */
#define ZCOM_PICK
#endif

#ifndef ZCOM_PICK
  #ifndef ZCOM_DEF
  #define ZCOM_DEF
  #endif
  #ifndef ZCOM_UTIL
  #define ZCOM_UTIL
  #endif
  #ifndef ZCOM_SS
  #define ZCOM_SS
  #endif
  #ifndef ZCOM_ENDN
  #define ZCOM_ENDN
  #endif
  #ifndef ZCOM_RNG
  #define ZCOM_RNG
  #endif
  #ifndef ZCOM_OPT
  #define ZCOM_OPT
  #endif
  #ifndef ZCOM_ARGOPT
  #define ZCOM_ARGOPT
  #endif
  #ifndef ZCOM_CFG
  #define ZCOM_CFG
  #endif
  #ifndef ZCOM_LOG
  #define ZCOM_LOG
  #endif
#endif

/* build dependencies */


#ifdef ZCOM_LOG
  #define ZCOM_UTIL
#endif

#ifdef ZCOM_CFG
  #define ZCOM_OPT
#endif

#ifdef ZCOM_ARGOPT
  #define ZCOM_OPT
#endif

#ifdef ZCOM_OPT
  #define ZCOM_SS
  #define ZCOM_UTIL
  #define ZCOM_DEF
#endif


/* manage storage class: static is the safer choice
   to avoid naming conclict.  Example:
   both m.c and n.c include this file,
   m.c --> m.o, n.c --> n.o, m.o+n.o --> a.out
   static is the only way to avoid naming conflict in this case.

   In case that this file is included multiple times,
   ZCOM_XFUNCS should be defined before the first inclusion,
   otherwise it won't be effective in deciding storage class. */
#ifndef STRCLS
  #ifndef ZCOM_XFUNCS
    #define STRCLS static
  #else
    #define STRCLS
  #endif
#endif

/* inline keyword */
#ifndef INLINE
  #if defined(__GNUC__) || defined(__xlC__)
    #define INLINE STRCLS __inline__
  #elif defined(__INTEL_COMPILER)
    #define INLINE STRCLS __inline
  #elif defined(_MSC_VER) || defined(__BORLANDC__)
    #define INLINE __inline STRCLS
  #elif defined(__STDC_VERSION__) && (STDC_VERSION__ >= 199901L)
    #define INLINE STRCLS inline
  #else
    #define INLINE STRCLS
  #endif
#endif

/* restrict keyword */
#ifndef RESTRICT
  #if (defined(__GNUC__) || defined(__INTEL_COMPILER) || defined(__xlC__))
    #define RESTRICT __restrict
  #elif defined(__STDC_VERSION__) && (STDC_VERSION__ >= 199901L)
    #define RESTRICT restrict
  #else
    #define RESTRICT
  #endif
#endif

/* macros with variable-length arguments */
#ifndef HAVEVAM
  #if (  (defined(__GNUC__) && (__GNUC__ >= 3))   \
      || (defined(__xlC__)  && (__xlC__ >= 0x0700)) \
      || (defined(_MSC_VER) && (_MSC_VER >= 1400)) )
    #define HAVEVAM 1
  #endif
#endif

#ifdef __INTEL_COMPILER
  #pragma warning(disable:981) /* unspecified order warning */
  #pragma warning(disable:177) /* unreferenced function */
  #pragma warning(disable:161) /* unrecognized #pragma, for omp */
#elif defined(__GNUC__) && (__GNUC__ >= 4 && __GNUC_MINOR__ >= 2)
  #pragma GCC diagnostic ignored "-Wunknown-pragmas"
  #pragma GCC diagnostic ignored "-Wvariadic-macros"
#endif

#ifdef _MSC_VER
  #pragma warning(disable:4127) /* conditional expression constant */
  #pragma warning(disable:4505) /* unreferenced function */
  #pragma warning(disable:4514) /* unreferenced inline */
  #pragma warning(disable:4710) /* not inlined */
  #include <stdio.h> /* suppress CRT _s functions warnings */
#endif

/* In addition to ZCOM_ABC, we have to define another macro ZCOM_ABC__
 * in order to avoid multiple inclusions.
 * A single ZCOM_ABC__ won't do because different module-set may be selected */

#ifdef  ZCOM_DEF
#ifndef ZCOM_DEF__
#define ZCOM_DEF__
#include <float.h>

/* define a real type */
#ifdef HAVE_REAL
  #ifndef HAVEREAL
  #define HAVEREAL HAVE_REAL
  #endif
#endif

#ifndef HAVEREAL
  #define HAVEREAL 1
  #define real double
#endif 

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif /* ZCOM_DEF__ */
#endif /* ZCOM_DEF */

#ifdef  ZCOM_UTIL
#ifndef ZCOM_UTIL__
#define ZCOM_UTIL__
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#ifndef xnew
#define xnew(x, n) \
  if (#n[0] != '1' && (n) <= 0) { \
    fprintf(stderr, "cannot allocate %d objects for %s\n", (int) (n), #x); \
    exit(1); \
  } else if ((x = calloc(n, sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %u\n", #x, (unsigned) (n)); \
    exit(1); }
#endif

#ifndef xrenew
#define xrenew(x, n) \
  if ((n) <= 0) { \
    fprintf(stderr, "cannot allocate %d objects for %s\n", (int) (n), #x); \
    exit(1); \
  } else if ((x = realloc(x, (n)*sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %u\n", #x, (unsigned) (n)); \
    exit(1); }
#endif

/* print an error message */
INLINE void perrmsg__(const char *file, int line, const char *why,
    const char *fmt, va_list args)
{
  fprintf(stderr, "error: ");
  vfprintf(stderr, fmt, args);
  if (fmt[strlen(fmt) - 1] != '\n')
    fprintf(stderr, "\n"); /* add a new line if needed */
  if (file != NULL) fprintf(stderr, "file: %s\n", file);
  if (line >= 0) fprintf(stderr, "line: %d\n", line);
  if (why != NULL && strcmp(why, "1") != 0)
    fprintf(stderr, "cond: %s\n", why);
}

#ifdef HAVEVAM

INLINE void perrmsg_(const char *file, int line, const char *why,
    int cond, const char *fmt, ...)
{
  if (cond) {
    va_list args;
    va_start(args, fmt);
    perrmsg__(file, line, why, fmt, args);
    va_end(args);
    exit(1);
  }
}

#define die_if(cond, fmt, ...) \
  perrmsg_(__FILE__, __LINE__, #cond, cond, fmt, ## __VA_ARGS__)
#define fatal(fmt, ...)  die_if(1, fmt, ## __VA_ARGS__)

#else /* !HAVEVAM */

#define PERRMSG__(c) {                        \
  if ((#c[0] == '1' && #c[1] == '\0') || c) { \
    va_list args;                             \
    va_start(args, fmt);                      \
    perrmsg__(NULL, -1, NULL, fmt, args);     \
    va_end(args);                             \
    exit(1);                                  \
  } }
INLINE void die_if(int cond, const char *fmt, ...) PERRMSG__(cond)
INLINE void fatal(const char *fmt, ...) PERRMSG__(1)
#undef PERRMSG__

#endif /* HAVEVAM */

#define xfopen(fp, fn, fmt, err) \
  if ((fp = fopen(fn, fmt)) == NULL) { \
    fprintf(stderr, "cannot open file %s\n", fn); err; }

INLINE int fexists(const char *fn)
{
  FILE *fp;
  if ((fp = fopen(fn, "r")) == NULL) return 0;
  else { fclose(fp); return 1; }
}

/* swap two variables */
#ifndef xtpswap
#define xtpswap(tp, x, y) { tp dum_; dum_ = (x); (x) = (y); (y) = dum_; }
#endif

#ifndef intswap
#define intswap(x, y) xtpswap(int, x, y)
#endif

#ifndef dblswap
#define dblswap(x, y) xtpswap(double, x, y)
#endif

INLINE int intmax(int x, int y) { return x > y ? x : y; }
INLINE int intmin(int x, int y) { return x < y ? x : y; }
/* confine x within [xmin, xmax] */
INLINE int intconfine(int x, int xmin, int xmax)
  { return x < xmin ? xmin : x > xmax ? xmax : x; }

INLINE int intsqr(int x) { return x * x; }

/* get the pair index from 0 to n*(n - 1)/2 - 1 */
INLINE int getpairindex(int i, int j, int n)
{
  die_if (i < 0 || i >= n || j < 0 || j >= n || i == j,
      "bad index error i %d, j %d, n %d\n", i, j, n);
  if (i > j) { int i1 = i; i = j; j = i1; }
  return n*i - (i + 1)*(i + 2)/2 + j;
}

/* return individual indices for a given pair index */
INLINE int parsepairindex(int id, int n, int *i, int *j)
{
  int i1, n1;
  die_if (id < 0 || id >= n*(n - 1)/2, "index %d too large for n %d\n", id, n);
  for (i1 = n - 2; i1 >= 0; i1--) {
    if (id >= (n1 = i1*n - i1*(i1 + 1)/2)) {
      *i = i1;
      *j = id - n1 + i1 + 1;
      return 0;
    }
  }
  return 1;
}

INLINE double dblmax(double x, double y) { return x > y ? x : y; }
INLINE double dblmin(double x, double y) { return x < y ? x : y; }
/* confine x within [xmin, xmax] */
INLINE double dblconfine(double x, double xmin, double xmax)
  { return x < xmin ? xmin : x > xmax ? xmax : x; }

INLINE double dblsqr(double x) { return x * x; }

/* sqrt(x*x + y*y) */
INLINE double dblhypot(double x, double y)
{
  double t;
  x = fabs(x);
  y = fabs(y);
  if (x <= 0.) return y;
  else if (y <= 0.) return x;
  if (x < y) t = x, x = y, y = t;
  t = y/x;
  return x*sqrt(1+t*t);
}

/* round x to a multiple dx  */
INLINE double dblround(double x, double dx)
{
  if (x*dx > 0) return dx * (int)(x/dx + (.5 - DBL_EPSILON));
  else return -dx * (int)(-x/dx + (.5 - DBL_EPSILON));
}

INLINE void dblcleararr(double *x, int n)
  { int i; for (i = 0; i < n; i++) x[i] = 0.0; }

#ifndef LNADD_DEFINED
#define LNADD_DEFINED
#define LN_BIG 50.0

/* log(exp(a) + exp(b)) */
INLINE double lnadd(double a, double b)
{
  double c;
  if (a < b) { c = a; a = b; b = c; } /* ensure a >= b */
  return ((c = a-b) > LN_BIG) ? a : a + log(1 + exp(-c));
}

/* log(exp(a) - exp(b)), only works for a > b */
INLINE double lndif(double a, double b)
{
  double c;
  die_if (a < b, "lndif: %g < %g\n", a, b);
  return ((c = a-b) > LN_BIG) ? a : a + log(1 - exp(-c));
}

/* log(exp(a)+b) */
INLINE double lnaddn(double a, double b)
{
  return (a > LN_BIG) ? a : a + log(1 + b*exp(-a));
}

#undef LN_BIG
#endif  /* LNADD_DEFINED */

#define cisalnum(c)   isalnum((unsigned char)(c))
#define cisalpha(c)   isalpha((unsigned char)(c))
#define cisdigit(c)   isdigit((unsigned char)(c))
#define cisxdigit(c)  isxdigit((unsigned char)(c))
#define cisprint(c)   isprint((unsigned char)(c))
#define cisspace(c)   isspace((unsigned char)(c))
#define cislower(c)   islower((unsigned char)(c))
#define cisupper(c)   isupper((unsigned char)(c))
#define ctolower(c)   (char) tolower((unsigned char)(c))
#define ctoupper(c)   (char) toupper((unsigned char)(c))

/* string manipulation */
#define ZSTR_XSPACEL  0x0001
#define ZSTR_XSPACER  0x0002
#define ZSTR_XSPACE   (ZSTR_XSPACEL|ZSTR_XSPACER)
#define ZSTR_COPY     0x0004
#define ZSTR_CAT      0x0008
#define ZSTR_CASE     0x0100
#define ZSTR_UPPER_   0x0200
#define ZSTR_UPPER    (ZSTR_CASE|ZSTR_UPPER_)
#define ZSTR_LOWER    ZSTR_CASE

/* remove leading and trailing spaces */
#define strip(s)  stripx(s, ZSTR_XSPACE)
#define lstrip(s) stripx(s, ZSTR_XSPACEL)
#define rstrip(s) stripx(s, ZSTR_XSPACER)
INLINE char *stripx(char *s, unsigned flags)
{
  char *p;

  if (flags & ZSTR_XSPACEL) { /* remove leading spaces */
    for (p = s; cisspace(*p); p++) ;
    if (*p == '\0') *s = '\0';
    else memmove(s, p, strlen(p)+1);
  }
  if (flags & ZSTR_XSPACER) /* remove trailing spaces */
    for (p = s + strlen(s) - 1; p >= s && cisspace(*p); p--)
      *p = '\0';
  return s;
}

/* in the follows, size_s means the buffer size of s, i.e., sizeof(s) for static strings */
/* copy the string and convert it to upper/lower case */
#define strcpy2u(s, t, size_s) strcnv(s, t, size_s - 1, ZSTR_COPY|ZSTR_UPPER)
#define strcpy2l(s, t, size_s) strcnv(s, t, size_s - 1, ZSTR_COPY|ZSTR_LOWER)
#define strcpy_sf(s, t, size_s) strcnv(s, t, size_s - 1, ZSTR_COPY)
#define substr(s, t, start, len) strcnv(s, t+start, len, ZSTR_COPY)
/* concatenate strings, the last parameter is the buffer size of s,
 * unlike strncat(), in which it's the number of characters from *t* to be copied.  */
#define strcat_sf(s, t, size_s) strcnv(s, t, size_s - 1, ZSTR_CAT)
/* safely copy/cat strings with case conversion
 * unlike strncpy(), s is always null-terminated on return: it copies at most
 * len nonblank characters, i.e., s[len] = '\0' for the longest output */
INLINE char *strcnv(char *s, const char *t, size_t len, unsigned flags)
{
  size_t i = 0, j;
  unsigned docase = flags & ZSTR_CASE, up = flags & ZSTR_UPPER_;

  if (len == 0 || s == NULL || t == NULL) return s;
  if (flags & ZSTR_CAT) while(s[i]) i++;
  for (j = 0; i < len; i++, j++) {
    if (docase && t[j]) {
      if (up) s[i] = (char) (unsigned char) toupper((unsigned char) t[j]);
      else    s[i] = (char) (unsigned char) tolower((unsigned char) t[j]);
    } else s[i] = t[j];
    if (t[j] == 0) break;
  }
  if (i == len) s[i] = '\0';
  if (flags & ZSTR_XSPACE) stripx(s, flags); /* call strip */
  return s;
}

/* compare strings without case */
#define strcmpnc(s, t) strncmpnc(s, t, -1)
INLINE int strncmpnc(const char *s, const char *t, int n)
{
  int i, cs, ct;

  if (s == NULL || t == NULL) return 0;
  for (i = 0; ; i++) {
    if (i >= n) return 0;
    cs = s[i];
    ct = t[i];
    if (cs == 0 || ct == 0) break;
    cs = toupper( (unsigned char) cs );
    ct = toupper( (unsigned char) ct );
    if (cs != ct) break;
  }
  return cs-ct;
}

#endif /* ZCOM_UTIL__ */
#endif /* ZCOM_UTIL */

#ifdef  ZCOM_SS
#ifndef ZCOM_SS__
#define ZCOM_SS__

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>

enum { SSCAT = 1, SSDELETE = 2, SSSHRINK = 3, SSSINGLE = 0x1000 };

#define ssnew(n)       sscpycatx(NULL, NULL, (n),    0)
#define ssdup(t)       sscpycatx(NULL, (t),   0,     0)
#define sscpy(s, t)    sscpycatx(&(s), (t),   0,     0)
#define sscat(s, t)    sscpycatx(&(s), (t),   0, SSCAT)
#define ssdel(s)       ssmanage((s), SSDELETE|SSSINGLE)
#define ssdelete(s)    { ssdel(s); (s)=NULL; }
#define ssshrink(s)    ssmanage((s), SSSHRINK|SSSINGLE)
#define ssdelall()     ssmanage(NULL, SSDELETE)
#define ssshrinkall()  ssmanage(NULL, SSHRINK)
#define ssfgets(s, pn, fp)    ssfgetx(&(s), (pn), '\n', (fp))
#define ssfgetall(s, pn, fp)  ssfgetx(&(s), (pn), EOF, (fp))

STRCLS int   ssmanage(char *, unsigned);
STRCLS char *sscpycatx(char **, const char *, size_t, unsigned);
STRCLS char *ssfgetx(char **, size_t *, int, FILE *fp);


#ifndef SSMINSIZ /* to override the block size, define it before inclusion */
#define SSMINSIZ 256 /* change this value to 1 for debugging */
#endif
#ifndef SSHASHBITS
#define SSHASHBITS 8
#endif
#define SSHASHSIZ   (1 << SSHASHBITS)
#define SSOVERALLOC 1
#define sscalcsize_(n) (((n)/SSMINSIZ + 1) * SSMINSIZ) /* size for n nonblank characters */

struct ssheader {
  size_t size;
  size_t hashval;
  struct ssheader *next;
};

static struct ssheader ssbase_[SSHASHSIZ] = {{ 0u, 0u, NULL }};

/* we use the string address instead of that of the pointer
 * to struct ssheader to compute the Hash value,
 * because the former is more frequently used in e.g. looking-up
 * */
static size_t sshashval_(const char *p)
{
  size_t val = (size_t) p * 1664525u + 1013904223u;
  return (val >> (sizeof(size_t)*8-SSHASHBITS)) & ((1<<SSHASHBITS)-1);
}

/*
 * return the *previous* header to the one that associates with s
 * first locate the list from the Hash value, then enumerate the linked list.
 * */
static struct ssheader *sslistfind_(const char *s)
{
  struct ssheader *hp, *head;

  if (s == NULL) return NULL;
  head = ssbase_ + sshashval_(s);
  if (head->next == NULL) return NULL; /* uninitialized head node */
  for (hp = head; hp->next != head; hp = hp->next)
    if ((char *)(hp->next + 1) == s)
      return hp;
  return NULL;
}

/*
 * simply add the entry h at the begining of the list
 * we do not accept a precalculated hash value,
 * since realloc might have changed it
 * */
static struct ssheader *sslistadd_(struct ssheader *h)
{
  struct ssheader *head;

  head = ssbase_ + sshashval_( (char *)(h+1) );
  if (head->next == NULL) /* initialize the base */
    head->next = head;
  h->next = head->next;
  head->next = h;
  return head;
}

/* remove hp->next */
static void sslistremove_(struct ssheader *hp, int f)
{
  struct ssheader *h = hp->next;

  hp->next = h->next;
  if (f) free(h);
}

/* (re)allocate memory for (*php)->next, update list, return the new string
 * n is the number of nonempty characters, obtained e.g. from strlen().
 * create a new header if *php is NULL, in this case, the first character
 * of the string is '\0'
 * */
static char *ssresize_(struct ssheader **php, size_t n, unsigned flags)
{
  struct ssheader *h = NULL, *hp;
  size_t size;

  if (php == NULL) {
    fprintf(stderr, "ssresize_: php is NULL, n = %u", (unsigned) n);
    return NULL;
  }

  /* we use the following if to assign hp and h, so the order is crucial */
  if ((hp = *php) == NULL || (h = hp->next)->size < n + 1 || !(flags & SSOVERALLOC)) {
    size = sscalcsize_(n);
    if (h == NULL || size != h->size) {
      /* since realloc will change the hash value of h
       * we have to remove the old entry first without free()
       * hp->next will be freed by realloc */
      if (hp != NULL)
        sslistremove_(hp, 0);
      if ((h = realloc(h, sizeof(*h)+size)) == NULL) {
        fprintf(stderr, "ssresize_: no memory for %u\n", (unsigned) size);
        return NULL;
      }
      if (hp == NULL) /* clear the first byte if we start from nothing */
        *(char *)(h + 1) = '\0';  /* h + 1 is the beginning of the string */
      *php = hp = sslistadd_(h);
      hp->next->size = size;
    }
  }
  return (char *)(hp->next + 1);
}

static void ssmanage_low_(struct ssheader *hp, unsigned opt)
{
  if (opt == SSDELETE)
    sslistremove_(hp, 1);
  else if (opt == SSSHRINK)
    ssresize_(&hp, strlen((char *)(hp->next+1)), 0);
}

/* delete a string, shrink memory, etc ... */
int ssmanage(char *s, unsigned flags)
{
  struct ssheader *hp, *head;
  unsigned opt = flags & 0xFF;
  size_t i;

  if (flags & SSSINGLE) {
    if (s == NULL || (hp = sslistfind_(s)) == NULL) {
      if (s) fprintf(stderr, "ssmanage: unknown address %p (%s)\n",  s, s);
      return -1;
    }
    ssmanage_low_(hp, opt);
  } else {
    for (i = 0; i < SSHASHSIZ; i++)
      for (hp = head = ssbase_+i; hp->next && hp->next != head; hp = hp->next)
        /* we must not operate on h itself, which renders the iterator h invalid */
        ssmanage_low_(hp, opt);
  }
  return 0;
}

/*
 * copy/cat t to *ps
 *
 * If (flags & SSCAT) == 0:
 * copy t to *ps, if ps is not NULL, and return the result
 * if ps or *ps is NULL, we return a string created from t
 *   *ps is set to the same value if ps is not NULL
 * otherwise, we update the record that corresponds to *ps
 *
 * minsize: to request a minimal size for the resulting buffer
 *
 * If flags & SSCAT:
 * append t after *ps. Equivalent to cpy if ps or *ps is NULL.
 * */
char *sscpycatx(char **ps, const char *t, size_t minsize, unsigned flags)
{
  struct ssheader *hp = NULL;
  size_t size = 0u, sizes = 0u;
  char *s = NULL, *p;

  /* both ps and *ps can be NULL, in which cases we leave hp as NULL */
  if (ps != NULL && (s = *ps) != NULL && (hp = sslistfind_(s)) == NULL) {
    fprintf(stderr, "sscpycatx: unknown address %p (%s)\n", s, s);
    return NULL;
  }
  if (t != NULL)
    while (t[size]) /* compute the length of t */
      size++;
  if (flags & SSCAT) {
    if (s != NULL)  /* s is also NULL, if ps itself is NULL */
      while (s[sizes]) /* compute the length of s */
        sizes++;
    size += sizes;
  }  /* sizes is always 0 in case of copying */
  if (size < minsize)
    size = minsize;
  if ((s = ssresize_(&hp, size, SSOVERALLOC)) == NULL) { /* change size */
    return NULL;
  }
  if (t != NULL)
    for (p = s + sizes; (*p++ = *t++); ) /* copy/cat the string */
      ;
  if (ps != NULL)
    *ps = s;
  return s;
}

/* get a string *ps from file fp
 * *ps can be NULL, in which case memory is allocated
 * *pn is number of characters read (including '\n', but not the terminal null)
 * delim is the '\n' for reading a singe line
 * */
char *ssfgetx(char **ps, size_t *pn, int delim, FILE *fp)
{
  size_t n, max;
  int c;
  char *s;
  struct ssheader *hp;

  if (ps == NULL || fp == NULL)
    return NULL;
  if ((s = *ps) == NULL) /* allocate an initial buffer if *ps is NULL */
    if ((s = sscpycatx(ps, NULL, 0, 0u)) == NULL)
      return NULL;
  if ((hp = sslistfind_(s)) == NULL) {
    fprintf(stderr, "ssfgetx: unknown address %p (%s)\n", s, s);
    return NULL;
  }
  max = hp->next->size-1;
  for (n = 0; (c = fgetc(fp)) != EOF; ) {
    if (n+1 > max) { /* request space for n+1 nonblank characters */
      if ((*ps = s = ssresize_(&hp, n+1, SSOVERALLOC)) == NULL)
        return NULL;
      max = hp->next->size - 1;
    }
    s[n++] = (char)(unsigned char) c;
    if (c == delim)
      break;
  }
  s[n] = '\0';
  if (pn != NULL)
    *pn = n;
  return (n > 0) ? s : NULL;
}
#endif /* ZCOM_SS__ */
#endif /* ZCOM_SS */

#ifdef  ZCOM_ENDN
#ifndef ZCOM_ENDN__
#define ZCOM_ENDN__


#include <stdio.h>
#include <string.h>

STRCLS int endn_system(void);
STRCLS size_t endn_fwrite(void *ptr, size_t size, size_t n, FILE *fp, int endn);
STRCLS size_t endn_fread(void *ptr, size_t size, size_t n, FILE *fp, int flip);
STRCLS int endn_rmatch(void *src, const void *ref, size_t size, FILE *fp);
STRCLS int endn_rmatchi(int *src, int iref, FILE *fp);


/* return the system endian, 1: big endian, 0: little endian */
int endn_system(void)
{
  unsigned feff = 0xFEFF; /* assume unsigned is at least 16-bit */
  unsigned char *p;

  p  = (unsigned char *) &feff;
  return (*p == 0xFF) ? 0 : 1;
}

/* change endianness in-place for n items of size in ptr */
INLINE void endn_flip(void *ptr, size_t size, size_t n)
{
  unsigned char *p = (unsigned char *) ptr, ch;
  size_t i, r, half = size/2;

  for (; n > 0; n--, p += size) {
    /* reverse bytes for each object */
    for (i = 0; i < half; i++) {
      r = size - i - 1;
      ch   = p[i];
      p[i] = p[r];
      p[r] = ch;
    }
  }
}

/* write data in ptr to file with a specific endian 'endn'
 * `ptr' is not const, because it needs to change its endian */
size_t endn_fwrite(void *ptr, size_t size, size_t n, FILE *fp, int endn)
{
  static int endsys = -1;

  /* initial determine the machine's endianess */
  if (endsys < 0) endsys = endn_system();
  if (endn == endsys) return fwrite(ptr, size, n, fp);

  endn_flip(ptr, size, n);
  n = fwrite(ptr, size, n, fp);
  endn_flip(ptr, size, n);
  return n;
}

/* read an object test object *src, compared with *ref
 * return 0 if they are identical without endian change
 * return 1 if changing the endianness of *src matches *ref
 * otherwise return -1 */
int endn_rmatch(void *src, const void *ref, size_t size, FILE *fp)
{
  if (1 != fread(src, size, 1, fp))
    return -1;
#ifdef ENDN_DBG
  if (size == sizeof(int))
    printf("A: 0x%X vs. 0x%X size = %u, cmp = %d\n",
      *(int *)src, *(int *)ref, (unsigned)size,
      memcmp(src, ref, size));
#endif
  if (memcmp(src, ref,  size) == 0)
    return 0;
  /* alter the endianness, and test again */
  endn_flip(src, size, 1);
#ifdef ENDN_DBG
  if (size == sizeof(int))
    printf("B: 0x%X vs. 0x%X size = %u, cmp = %d\n",
      *(int *)src, *(int *)ref, (unsigned)size,
      memcmp(src, ref, size));
#endif
  return (memcmp(src, ref, size) == 0) ? 1 : -1;
}

/* special case of endn_rmatchi for integer, convenient because
 * iref could be e.g. sizeof(int), which has no address */
int endn_rmatchi(int *src, int iref, FILE *fp)
{
  return endn_rmatch(src, &iref, sizeof(int), fp);
}

/* read data from file to ptr with endianness changed if 'flip' is 1
 * flip can be initialized by calling endn_rmatch() for a test object */
size_t endn_fread(void *ptr, size_t size, size_t n, FILE *fp, int flip)
{
  n = fread(ptr, size, n, fp);
  if (flip) endn_flip(ptr, size, n);
  return n;
}

#endif /* ZCOM_ENDN__ */
#endif /* ZCOM_ENDN */


#ifdef  ZCOM_RNG
#ifndef ZCOM_RNG__
#define ZCOM_RNG__

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
  #include <inttypes.h>
#elif defined(_MSC_VER) || defined(__BORLANDC__)
  typedef unsigned uint32_t;
  typedef unsigned __int64 uint64_t;
#else
  #include <inttypes.h>
#endif

#ifndef PRIu32
  #if (defined(_MSC_VER) && (_MSC_VER >= 1300)) || defined(__BORLANDC__)
    #define PRIu32 "I32u"
  #else
    #define PRIu32 "u"
  #endif
#endif

#ifndef PRIu64
  #if defined(_MSC_VER) || defined(__BORLANDC__)
    #define PRIu64 "I64u"
  #else
    #define PRIu64 "llu"
  #endif
#endif

#define rand32()  mtrand()
#define rnd0()    ((1.0/4294967296.0) * rand32()) /* double, [0, 1) */

#define MTFILE    "MTSEED"  /* default file */
#define MTSEED    time(0)   //5489UL    /* default seed */
STRCLS int mtsave(const char *fname);
STRCLS int mtload(const char *fname, uint32_t seed);
STRCLS uint32_t mtrand(void);
STRCLS double grand0(void);

/* metropolis acceptance probability rnd0() < exp(r), assuming r > 0 */
INLINE int metroacc0(double r) { r = exp(r); return rnd0() < r; }

/* metropolis acceptance probability rnd0() < exp(bet * r), assuming bet > 0
 * defined as a macro, in case r is an integer */
#define metroacc(r, bet) ((r >= 0) ? 1 : metroacc0(r * bet))


/* Mersenne Twister was developped by Makoto Matsumoto and Takuji Nishimura */
#define MT_N 624
#define MT_M 397
#define MT_UMASK 0x80000000UL /* most significant w-r bits */
#define MT_LMASK 0x7fffffffUL /* least significant r bits */

STRCLS int mtidx_ = -1; /* index in mt_, -1: uninitialized */
uint32_t mt_[MT_N]; /* array for the mt state vector */

/* save the current mt state to file */
int mtsave(const char *fname)
{
  FILE *fp;
  int k;

  if (mtidx_ < 0) return 1; /* RNG was never used, so it cannot be saved */
  if (fname == NULL) fname = MTFILE;
  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot save to %s.\n", fname);
    return 1;
  }
  fprintf(fp, "MTSEED\n%d\n", mtidx_);
  for (k = 0; k < MT_N; k++) fprintf(fp, "%"PRIu32"\n", mt_[k]);
  fclose(fp);
  return 0;
}

/* load mt state from `fname', or if it fails, use `seed' to initialize mt  */
int mtload(const char *fname, uint32_t seed)
{
  static char s[64];
  int k, z, err = 1;
  FILE *fp;

  if (fname == NULL) fname = MTFILE;
  if ((fp = fopen(fname, "r")) != NULL) { /* try to load from file */
    if (fgets(s, sizeof s, fp) == NULL) {
      fprintf(stderr, "%s is empty\n", fname);
    } else if (strncmp(s, "MTSEED", 6) != 0) { /* to check the first line */
      fprintf(stderr, "mtrand: corrupted file.\n");
    } else if (fscanf(fp, "%d", &mtidx_) != 1) {
      fprintf(stderr, "no index in %s\n", fname);
    } else {
      if (mtidx_ < 0) mtidx_ = MT_N; /* request updating */
      for (z = 1, k = 0; k < MT_N; k++) {
        if (fscanf(fp, "%"PRIu32, &mt_[k]) != 1) break;
        if (mt_[k] != 0) z = 0; /* a non-zero number */
      }
      if (k != MT_N) fprintf(stderr, "%s incomplete %d/%d\n", fname, k, MT_N);
      else err = z; /* clear error, if array is nonzero */
    }
    fclose(fp);
  }

  if (err) { /* initialize from seed */
    if (seed == 0) seed = MTSEED;
		fprintf(stderr, "mtload: cannot init from file:\"%s\". Will init from seed: %d directly.\n", fname, seed);
		mt_[0] = seed & 0xffffffffUL;
    for (k = 1; k < MT_N; k++) /* the final mask is for 64-bit machines */
      mt_[k] = (1812433253UL * (mt_[k-1] ^ (mt_[k-1]>>30)) + k) & 0xffffffffUL;
    mtidx_ = MT_N; /* request updating */
  }
	else
		fprintf(stderr, "Successfully loaded\"%s\"!\n", fname);
  
	return (mtidx_ < 0);
}

/* return an unsigned random number */
uint32_t mtrand(void)
{
  uint32_t x;
  static const uint32_t mag01[2] = {0, 0x9908b0dfUL}; /* MATRIX_A */
  int k;

  if (mtidx_ < 0) mtload(NULL, 0);
  if (mtidx_ >= MT_N) { /* generate MT_N words at one time */
    for (k = 0; k < MT_N - MT_M; k++) {
      x = (mt_[k] & MT_UMASK) | (mt_[k+1] & MT_LMASK);
      mt_[k] = mt_[k+MT_M] ^ (x>>1) ^ mag01[x&1UL];
    }
    for (; k < MT_N-1; k++) {
      x = (mt_[k] & MT_UMASK) | (mt_[k+1] & MT_LMASK);
      mt_[k] = mt_[k+(MT_M-MT_N)] ^ (x>>1) ^ mag01[x&1UL];
    }
    x = (mt_[MT_N-1] & MT_UMASK) | (mt_[0] & MT_LMASK);
    mt_[MT_N-1] = mt_[MT_M-1] ^ (x>>1) ^ mag01[x&1UL];
    mtidx_ = 0;
  }
  x = mt_[ mtidx_++ ];
  /* tempering */
  x ^= (x >> 11);
  x ^= (x <<  7) & 0x9d2c5680UL;
  x ^= (x << 15) & 0xefc60000UL;
  x ^= (x >> 18);
  return x;
}

#undef MT_N
#undef MT_M
#undef MT_UMASK
#undef MT_LMASK

/* Gaussian distribution with zero mean and unit variance
 * using ratio method */
double grand0(void)
{
  double x, y, u, v, q;
  do {
    u = 1 - rnd0();
    v = 1.7156*(rnd0() - .5);  /* >= 2*sqrt(2/e) */
    x = u - 0.449871;
    y = fabs(v) + 0.386595;
    q = x*x  + y*(0.196*y - 0.25472*x);
    if (q < 0.27597) break;
  } while (q > 0.27846 || v*v > -4*u*u*log(u));
  return v/u;
}

#endif /* ZCOM_RNG__ */
#endif /* ZCOM_RNG */


#ifdef  ZCOM_OPT
#ifndef ZCOM_OPT__
#define ZCOM_OPT__
#include <stdio.h>

/* option either from arguments or configuration */
typedef struct {
  int isopt; /* is option or argument */
  char ch; /* single letter option flag */
  const char *sflag; /* long string flag */
  const char *key; /* key */

  const char *val; /* raw string from command line */
  const char *desc; /* description */
  const char *fmt; /* sscanf format */
  const char *pfmt; /* printf format, NULL: to guess */
  void *ptr; /* address to the target variable */
  unsigned flags;
} opt_t;

#define OPT_MUST     0x0001  /* a mandatory argument or option */
#define OPT_SWITCH   0x0002  /* an option is a switch */
#define OPT_SET      0x0004  /* an argument/option is set */

/* translate string values to actual ones through sscanf() */
INLINE int opt_getval(opt_t *o)
{
  const char *fmt = o->fmt;

  if (fmt == NULL || fmt[0] == '\0') { /* raw string assignment */
    *(const char **)o->ptr = o->val;
  } else if (strcmp(fmt, "%s") == 0) {
    sscpy( *(char **)o->ptr, o->val);
  } else { /* call sscanf */
    if (strcmp(fmt, "%r") == 0) /* real */
      fmt = (sizeof(real) == sizeof(float)) ? "%f" : "%lf";
    if (1 != sscanf(o->val, fmt, o->ptr)) {
      fprintf(stderr, "Error: unable to convert a value for [%s] as fmt [%s], raw string: [%s]\n",
          o->desc, fmt, o->val);
      return 1;
    }
  }
  return 0;
}

/* set properties of an option: fmt = "%b" for a switch */
INLINE void opt_set(opt_t *o, const char *sflag, const char *key,
    const char *fmt, void *ptr, const char *desc)
{
  o->ch = '\0';
  if (key) {
    o->isopt = 2;
  } else if (sflag) { /* option */
    o->isopt = 1;
    o->ch = (char) ( sflag[2] ? '\0' : sflag[1] ); /* no ch for a long flag */
  } else { /* argument */
    o->isopt = 0;
  }
  o->sflag = sflag;
  o->key = key;
  o->flags = 0;
  die_if (ptr == NULL, "null pass to argopt with %s: %s\n", sflag, desc);
  o->ptr = ptr;
  if (fmt == NULL) fmt = "";
  if (fmt[0] == '!') {
    fmt++;
    o->flags |= OPT_MUST;
  }
  if (strcmp(fmt, "%b") == 0) {
    fmt = "%d";
    o->flags |= OPT_SWITCH;
  }
  o->fmt = fmt;
  o->pfmt = NULL;
  o->desc = desc;
}

/* print the value of o->ptr */
INLINE void opt_printptr(opt_t *o)
{
  const char *fmt;

  for (fmt = o->fmt; *fmt && *fmt != '%'; fmt++) ;
#define ELIF_PF_(fm, fmp, type) else if (strcmp(fmt, fm) == 0) printf((o->pfmt ? o->pfmt : fmp), *(type *)o->ptr)
  if (fmt == NULL || *fmt == '\0') printf("%s", (*(char **)o->ptr) ? (*(char **)o->ptr) : "NULL");
  ELIF_PF_("%b", "%d", int); /* switch */
  ELIF_PF_("%d", "%d", int);
  ELIF_PF_("%u", "%u", unsigned);
  ELIF_PF_("%x", "0x%x", unsigned);
  ELIF_PF_("%ld", "%ld", long);
  ELIF_PF_("%lu", "%lu", unsigned long);
  ELIF_PF_("%lx", "0x%lx", unsigned long);
#if 0  /* C99 only */
  ELIF_PF_("%lld", "%lld", long long);
  ELIF_PF_("%llu", "%llu", unsigned long long);
  ELIF_PF_("%llx", "0x%llx", unsigned long long);
#endif
  ELIF_PF_("%f", "%g", float);
  ELIF_PF_("%lf", "%g", double);
  ELIF_PF_("%r", "%g", real);
  else printf("unknown %s-->%%d: %d\n", fmt, *(int *)o->ptr);
#undef ELIF_PF_
}

/* search an option list, return an option whose variable address is p */
INLINE opt_t *opt_find(opt_t *ls, int n, const void *p)
{
   int i;
   for (i = 0; i < n; i++) if (ls[i].ptr == p) return ls + i;
   return NULL;
}

/* search an option list to see if an option is explicitly set */
INLINE int opt_isset(opt_t *ls, int n, const void *p, const char *var)
{
  opt_t *o = opt_find(ls, n, p);
  die_if (!o, "cannot find var %s, ptr %p\n", var, p);
  return o->flags & OPT_SET ? 1 : 0;
}
#endif /* ZCOM_OPT__ */
#endif /* ZCOM_OPT */

#ifdef  ZCOM_ARGOPT
#ifndef ZCOM_ARGOPT__
#define ZCOM_ARGOPT__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
  int nopt;
  opt_t *opts;
  const char *prog;
  const char *desc;
  const char *author;
  const struct tm *tm; /* compilation time */
  int version;
  unsigned flags;
  int dum_[4]; /* space holder */
} argopt_t;

#define ARGOPT_MUST     OPT_MUST    /* mandatory argument or option, format starts with ! */
#define ARGOPT_SWITCH   OPT_SWITCH  /* format "%b" */
#define ARGOPT_SET      OPT_SET
#define ARGOPT_LONGOPT  0x0010  /* always assume long format, e.g., -maxh */

STRCLS argopt_t *argopt_open(unsigned flags);
STRCLS void argopt_close(argopt_t *ao);
#define argopt_regarg(ao, fmt, ptr, desc) argopt_add(ao, NULL, fmt, ptr, desc)
#define argopt_regopt argopt_add
#define argopt_reghelp argopt_addhelp
#define argopt_regversion argopt_addversion
STRCLS int argopt_add(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc);
STRCLS void argopt_parse(argopt_t *ao, int argc, char **argv);

#define argopt_addhelp(ao, sflag) argopt_add(ao, sflag, "%b", ao->dum_, "$HELP")
#define argopt_addversion(ao, sflag) argopt_add(ao, sflag, "%b", ao->dum_, "$VERSION")

#define argopt_getopt(ao, p) opt_find(ao->opts, ao->nopt, p)
#define argopt_getarg argopt_getopt

/* test if argument/option is explicitly set */
#define argopt_set(ao, var) opt_isset(ao->opts, ao->nopt, &var, #var)


/* initialize the argument structure */
argopt_t *argopt_open(unsigned flags)
{
  argopt_t *ao;
  time_t tmcmpl;

  xnew(ao, 1);
  ao->flags = flags;
  ao->nopt = 0;
  ao->opts = NULL;
  tmcmpl = time(NULL);
  ao->tm = localtime( &tmcmpl );
  memset(ao->dum_, '\0', sizeof(ao->dum_));
  return ao;
}

void argopt_close(argopt_t *ao)
{
  if (ao->opts) { free(ao->opts); ao->opts = NULL; }
  free(ao);
}

/* print version and die */
static void argopt_version(argopt_t *ao)
{
  printf("%s: %s, version %d\n",
      ao->prog, ao->desc ? ao->desc : "", ao->version);
  if (ao->author && ao->tm)
    printf("Copyright (c) %s %d\n", ao->author, ao->tm->tm_year + 1900);
  argopt_close(ao);
  exit(1);
}

/* print help message and die */
static void argopt_help(argopt_t *ao)
{
  int i, len, maxlen;
  opt_t *o;
  const char *sysopt[2] = {"print help message", "print version"}, *desc;

  printf("%s, version %d",
      ao->desc ? ao->desc : ao->prog, ao->version);
  if (ao->author && ao->tm)
    printf(", Copyright (c) %s %d", ao->author, ao->tm->tm_year + 1900);
  printf("\nUSAGE\n  %s [OPTIONS]", ao->prog);
  for (i = 0; i < ao->nopt; i++) {
    const char *bra = "", *ket = "";
    o = ao->opts + i;
    if (o->isopt) continue;
    if (o->flags & OPT_MUST) {
      if (strchr(o->desc, ' '))
        bra = "{", ket = "}";
    } else
      bra = "[", ket = "]";
    printf(" %s%s%s", bra, o->desc, ket);
  }
  printf("\n");

  printf("OPTIONS:\n") ;
  for (maxlen = 0, i = 0; i < ao->nopt; i++) { /* compute the longest option */
    if (!ao->opts[i].isopt) continue;
    len = strlen(ao->opts[i].sflag);
    if (len > maxlen) maxlen = len;
  }
  for (i = 0; i < ao->nopt; i++) {
    o = ao->opts + i;
    if (!o->isopt) continue;
    desc = o->desc;
    if (strcmp(desc, "$HELP") == 0)
      desc = sysopt[0];
    else if (strcmp(desc, "$VERSION") == 0)
      desc = sysopt[1];
    printf("  %-*s : %s%s", maxlen, o->sflag,
        (!(o->flags & OPT_SWITCH) ? "followed by " : ""), desc);
    if (o->ptr && o->ptr != ao->dum_) { /* print default values */
      printf(", default: ");
      opt_printptr(o);
    }
    printf("\n");
  }
  argopt_close(ao);
  exit(1);
}

/* register an argument or option
 * sflag: string flag, or NULL for an argument
 * fmt: sscanf() format string, "%b" for a switch, "%r" for real
 * return the index */
int argopt_add(argopt_t *ao, const char *sflag,
    const char *fmt, void *ptr, const char *desc)
{
  opt_t *o;
  int n;

  n = ao->nopt++;
  xrenew(ao->opts, ao->nopt);
  o = ao->opts + n;
  opt_set(o, sflag, NULL, fmt, ptr, desc);
  return n;
}

/* main parser of arguments */
void argopt_parse(argopt_t *ao, int argc, char **argv)
{
  int i, j, k, ch, acnt = 0;
  opt_t *ol = ao->opts;

  ao->prog = argv[0];
  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') { /* it's an argument */
      while (ol[acnt].isopt && acnt < ao->nopt) acnt++;
      if (acnt >= ao->nopt) argopt_help(ao);
      ol[acnt].val = argv[i];
      ol[acnt].flags |= OPT_SET;
      if (0 != opt_getval(ol + acnt))
        argopt_help(ao);
      ++acnt;
      continue;
    }

    /* it's an option, loop for abbreviated form "-abc" == "-a -b -c" */
    for (j = 1; (ch = argv[i][j]) != '\0'; j++) {
      int islong = (j == 1 && argv[i][1] == '-') | (ao->flags & ARGOPT_LONGOPT);

      if (islong) { /* compare against long options */
        for (k = 0; k < ao->nopt; k++)
          if (ol[k].isopt &&
              strncmp(argv[i], ol[k].sflag, strlen(ol[k].sflag)) == 0)
            break;
      } else { /* compare against short options */
        for (k = 0; k < ao->nopt; k++)
          if (ol[k].isopt && ch == ol[k].ch)
            break;
      }
      if (k >= ao->nopt) {
        fprintf(stderr, "cannot handle option [%s]\n", argv[i]);
        argopt_help(ao);
      }

      if (ol[k].desc[0] == '$') { /* system commands */
        if (strcmp(ol[k].desc, "$HELP") == 0)
          argopt_help(ao);
        else if (strcmp(ol[k].desc, "$VERSION") == 0)
          argopt_version(ao);
      }

      if (ol[k].flags & OPT_SWITCH) {
        ol[k].flags |= OPT_SET;
        *(int *)ol[k].ptr = 1;
        if (islong) break; /* go to the next argument argv[i+1] */
      } else { /* look for the additional argument for this */
        int hasv = 0;
        if (islong) { /* e.g., --version=11 */
          j = strlen(ol[k].sflag);
          if (argv[i][ j ] == '=') {
            ol[k].val = argv[i] + j + 1;
            hasv = 1;
          }
        } else { /* e.g., -n8 */
          if (argv[i][++j]) {
            ol[k].val = argv[i] + j;
            hasv = 1;
          }
        }

        if (!hasv) { /* --version 11 or -n 8 */
          if (++i >= argc) {
            printf("%s(%s) requires an argument!\n", ol[k].sflag, argv[i-1]);
            argopt_help(ao);
          }
          ol[k].val = argv[i];
        }
        ol[k].flags |= OPT_SET;
        if (0 != opt_getval(ol+k)) argopt_help(ao);
        break; /* go to the next argument argv[i+1] */
      }
    } /* end of short option loop */
  }
  /* check if we have all mandatory arguments and options */
  for (i = 0; i < ao->nopt; i++) {
    if ((ol[i].flags & OPT_MUST) && !(ol[i].flags & OPT_SET)) {
      printf("Error: missing %s %s: %s\n\n",
          ol[i].isopt ? "option" : "argument", ol[i].sflag, ol[i].desc);
      argopt_help(ao);
    }
  }
}

#endif /* ZCOM_ARGOPT__ */
#endif /* ZCOM_ARGOPT */

#ifdef  ZCOM_CFG
#ifndef ZCOM_CFG__
#define ZCOM_CFG__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct {
  char *key, *val;
  int used;
} cfgent_t; /* line from cfg file */

typedef struct {
  char *buf;      /* the entire configuration file */
  int nent;       /* number of entries */
  cfgent_t *ents; /* entries */
  int nopt;       /* number of user-requested options */
  opt_t *opts;    /* user-requested options */
} cfg_t;
typedef cfg_t cfgdata_t;

#define CFG_CHECKUSE 0x0100
#define CFG_VERBOSE  0x1000

STRCLS cfg_t *cfg_open(const char *fn);
STRCLS void cfg_close(cfg_t *cfg);
STRCLS int cfg_add(cfg_t *cfg, const char *key, const char *fmt, void *ptr, const char *desc);
STRCLS int cfg_match(cfg_t *cfg, unsigned flags);

#define cfg_set(cfg, var) opt_isset(cfg->opts, cfg->nopt, &var, #var)

/* old style functions */
#define cfgopen(fn) cfg_open(fn)
#define cfgclose(cfg) cfg_close(cfg)
/* Read the value of a given variable from the current configuration file,
 * the name of variable is given by `key',
 * If the key is matched, its value is saved to `*var' through sscanf,
 *   otherwise, the content in *var is not modified.
 * If the function succeeds, it returns 0.
 * In case fmt is "%s", (*var) is a string, or a pointer to char.
 *   The space for (*var) will be managed through sscpy. */
INLINE int cfgget(cfg_t *cfg, void *var, const char *key, const char *fmt)
{
  int i;

  for (i = 0; i < cfg->nent; i++) {
    cfgent_t *ent = cfg->ents + i;
    if (ent->key != NULL && strcmp(ent->key, key) == 0) {
      if (strcmp(fmt, "%s") == 0) { /* string */
        sscpy( *(char **)var, ent->val); /* make a copy and return */
        return 0;
      } else /* use sscanf for other cases, like int, float,... */
        return EOF == sscanf(ent->val, fmt, var) ? 2 : 0;
    }
  }
  return 1; /* no match */
}

/* load the whole configuration file into memory, parse it to entries */
cfg_t *cfg_open(const char *fn)
{
  cfg_t *cfg;
  cfgent_t *ent;
  FILE *fp;
  size_t i, j, n, size = 0;
  char *p, *q;

  xnew(cfg, 1);

  xfopen(fp, fn, "r", return NULL);
  if (ssfgetall(cfg->buf, &size, fp) == NULL) {
    fprintf(stderr, "error reading file %s\n", fn);
    return NULL;
  }
  sscat(cfg->buf, "\n"); /* in case the file is not ended by a new line, we add one */
  fclose(fp);

  /* count the number of lines (before allocating the key-table) */
  for (p = cfg->buf, i = 0, n = 0; i < size; i++) {
    if (p[i] == '\n' || p[i] == '\r') {
      if (i > 0 && p[i-1] == '\\') {
        /* multiple-line splicing by replacing cr, lf with spaces
           parse should be aware of these additional spaces */
        p[i-1] = p[i] = ' ';
      } else {
        p[i] = '\0';
        n++;
      }

      /* replace following CR LF by spaces for efficiency
         as the size of the key table == the number of blank lines */
      for (j = i+1; j < size && cisspace(p[j]); j++) p[j] = ' ';
    }
  }

  xnew(cfg->ents, n);

  /* load lines into the keytable */
  for (p = cfg->buf, j = 0, i = 0; i < size; i++) {
    if (cfg->buf[i] == '\0') {
      cfg->ents[j++].key = p;
      if (j >= n) break;
      p = cfg->buf + i + 1;
    }
  }
  n = j;

  /* parse each line to a key-value pair */
  for (j = 0; j < n; j++) {
    ent = cfg->ents + j;
    p = ent->key;
    strip(p);

    /* skip a blank or comment line */
    if (p[0] == '\0' || strchr("#%!;", p[0]) != NULL) {
      ent->key = NULL;
      continue;
    }

    /* remove trailing space and ';' */
    for (q = p + strlen(p) - 1; q >= p && (cisspace(*q) || *q == ';'); q--)
      *q = '\0';

    if ((q = strchr(p, '=')) == NULL) { /* skip a line without '=' */
      ent->key = NULL;
      continue;
    }
    *q = '\0';
    ent->val = q + 1;
    strip(ent->key);
    strip(ent->val);
  }
  cfg->nent = (int) n;
  cfg->nopt = 0;
  xnew(cfg->opts, 1);
  return cfg;
}

void cfg_close(cfg_t *cfg)
{
  ssdelete(cfg->buf);
  free(cfg->ents);
  free(cfg->opts);
  memset(cfg, 0, sizeof(*cfg));
  free(cfg);
}

/* register an option request, return the index */
int cfg_add(cfg_t *cfg, const char *key, const char *fmt, void *ptr, const char *desc)
{
  int n = cfg->nopt++;
  opt_t *o;
  xrenew(cfg->opts, cfg->nopt);
  o = cfg->opts + n;
  opt_set(o, NULL, key, fmt, ptr, desc);
  return n;
}

/* match requested options with entries in cfg file
 * returns 0 if successful */
int cfg_match(cfg_t *cfg, unsigned flags)
{
  int i, j, ret = 0, verbose = flags & CFG_VERBOSE, must = flags & OPT_MUST;
  opt_t *o;
  cfgent_t *ent;

  for (i = 0; i < cfg->nopt; i++) {
    o = cfg->opts + i;
    for (j = 0; j < cfg->nent; j++) {
      ent = cfg->ents + j;
      if (ent->key != NULL && strcmp(ent->key, o->key) == 0) {
        ent->used = 1;
        o->flags |= OPT_SET;
        o->val = ent->val;
        opt_getval(o);
        break;
      }
    }
    if (!(o->flags & OPT_SET) && (must || verbose)) {
      printf("cfg: %s not set, default: ", o->key);
      opt_printptr(o);
      printf("\n");
      if (must) ret = 1;
    }
  }

  if (flags & CFG_CHECKUSE) {
    for (j = 0; j < cfg->nent; j++) {
      ent = cfg->ents + j;
      if (ent->key != NULL && !ent->used && verbose)
        printf("cfg: unused entry: %s = %s\n", ent->key, ent->val);
    }
  }
  return ret;
}

#endif /* ZCOM_CFG__ */
#endif /* ZCOM_CFG */

#ifdef  ZCOM_LOG
#ifndef ZCOM_LOG__
#define ZCOM_LOG__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

typedef struct {
  FILE *fp;
  const char *fname;
  unsigned flags;
} logfile_t;

#define LOG_WRITESCREEN  0x01
#define LOG_FLUSHAFTER   0x02
#define LOG_NOWRITEFILE  0x10
#define LOG_APPEND       0x80

INLINE logfile_t *log_open(const char *fn)
{
  logfile_t *log;

  xnew(log, 1);
  if (fn == NULL) fn = "LOG";
  log->fname = fn;
  log->flags = 0;
	log->flags |= LOG_APPEND; 
  return log;
}

INLINE int log_printf(logfile_t *log, char *fmt, ...)
{
  va_list args;

  if (log == NULL) return 1;
  if (log->fp == NULL) {
    const char *aw = (log->flags & LOG_APPEND) ? "a" : "w";
    xfopen(log->fp, log->fname, aw, return 1);
  }
  if ((log->flags & LOG_NOWRITEFILE) == 0) {
    va_start(args, fmt);
    vfprintf(log->fp, fmt, args);
    va_end(args);
  }
  if (log->flags & LOG_WRITESCREEN) {
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }
  if (log->flags & LOG_FLUSHAFTER)
    fflush(log->fp);
  return 0;
}

INLINE void log_close(logfile_t *log)
{
  if (log == NULL) return;
  if (log->fp != NULL) { fclose(log->fp); log->fp = NULL; }
  free(log);
}


/* close & reopen log file to make sure that stuff is written to disk */
INLINE int log_hardflush(logfile_t *log)
{
  if (log->fp == NULL || log->fname == NULL) return 1;
  fclose(log->fp);
  xfopen(log->fp, log->fname, "a", return 1);
  return 0;
}

#if defined(HAVEVAM) && defined(NEED_WTRACE)
STRCLS logfile_t log_stock_[1] = {{ NULL, "TRACE", 0 }};
#define wtrace(fmt, ...) { \
  if (fmt) log_printf(log_stock_, fmt, ##__VA_ARGS__); \
  else if (log_stock_->fp) { fclose(log_stock_->fp); log_stock_->fname = NULL; } }
#endif

#endif /* ZCOM_LOG__ */
#endif /* ZCOM_LOG */


