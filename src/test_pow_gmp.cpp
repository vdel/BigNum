#include "monitor.h"

#include <stdlib.h>
#include <stdio.h>
#include "iddlib.h"
#include "gmp.h"
#include "timing.h"

extern long call_constructor, call_clean_Ht;
extern long call_T, call_long_T;
extern long call_COMP, call_C;
extern long call_OF, saved_call_OF;
extern long call_DEC, call_INC, call_long_INC;
extern long call_X, saved_call_X;
extern long call_RMSB, call_AMSB;
extern long call_TWICE, saved_call_TWICE;
extern long call_A, saved_call_A, call_long_A;
extern long call_P, saved_call_P, call_long_P;
extern long call_NEG;
extern long call_LSH, call_RSH;
extern long call_POW, call_Y, call_HALF;


int cputime (void);

int
main (int argc, char *argv[])
{
  gmp_randstate_t rs;
  mpz_t x, y, z, res;
  mpz_ptr xptr, yptr;
  unsigned long int m, n, i, niter, t0, ti;
  double t;   
  uIDD x2, y2, z2;

  gmp_randinit_default (rs);

  mpz_init (x);
  mpz_init (y);
  mpz_init (z);
  mpz_init (res);  

  if (argc == 2)
    {
      m = atoi (argv[1]);
      n = m;
      mpz_urandomb (x, rs, m);
      xptr = x;
      yptr = x;
      x2 = x;
      y2 = x;
    }
  else if (argc == 3)
    {
      m = atoi (argv[1]);
      n = atoi (argv[2]);
      mpz_urandomb (x, rs, m);
      mpz_urandomb (y, rs, n);
      xptr = x;
      yptr = y;
      x2 = x;
      y2 = y;
    }
  else
    {
      fprintf (stderr, "usage: %s m n\n", argv[0]);
      fprintf (stderr, "  where m and n are number of bits in numbers tested\n");
      return -1;
    }

  double vm1, res1;
  process_mem_usage(vm1,res1);
  
  // Power
  double pow_per_sec_gmp;
  printf ("Calibrating CPU speed...");  fflush (stdout);
  TIME (t, mpz_pow_ui (z, xptr, 100));
  printf ("done\n");
  niter = 1 + (unsigned long) (1e4 / t);
  
  printf ("Computing the 100-power of a %lu-bit number %lu times with IDD...", m, niter);
  fflush (stdout);     
  t0 = cputime ();
  for (i = niter; i > 0; i--)
    mpz_pow_ui (z, xptr, 100);
  ti = cputime () - t0;
  printf ("done!\n");
  pow_per_sec_gmp = 1000.0 * niter / ti;

  double vm2, res2;
  process_mem_usage(vm2,res2);

  printf ("RESULTS:\n");
  printf("GMP: %f 100-power per second\n", pow_per_sec_gmp);
  printf("Memory used: %dKo\n\n", (int)(res2-res1));
  return 0;
}

/* Return user CPU time measured in milliseconds.  */
#if !defined (__sun) \
    && (defined (USG) || defined (__SVR4) || defined (_UNICOS) \
	|| defined (__hpux))
#include <time.h>

int
cputime ()
{
  return (int) ((double) clock () * 1000 / CLOCKS_PER_SEC);
}
#else
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>

int
cputime ()
{
  struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}
#endif
