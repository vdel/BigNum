#include "monitor.h"

#include <stdlib.h>
#include <stdio.h>
#include "iddlib.h"
#include "gmp.h"
#include "timing.h"

#define HIDE_RESULTS

extern int call_constructor, call_clean_Ht;
extern int call_T, call_int_T;
extern int call_COMP, call_C;
extern int call_OF, saved_call_OF;
extern int call_DEC, call_INC, call_int_INC;
extern int call_X, saved_call_X;
extern int call_RMSB, call_AMSB;
extern int call_TWICE, saved_call_TWICE;
extern int call_A, saved_call_A, call_int_A;
extern int call_P, saved_call_P, call_int_P;
extern int call_NEG;
extern int call_LSH, call_RSH;
extern int call_POW, call_Y, call_HALF;


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
      mpz_set(y, x);
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
 
  double mul_per_sec_idd;   
  printf ("Calibrating CPU speed...");  fflush (stdout);
  TIME (t, z2 = x2*y2);
  printf ("done\n");
  niter = 1 + (unsigned long) (1e4 / t);      
  
  IDDContent::clean_Ht();
  if (argc == 2)
    printf ("Squaring a %lu-bit number %lu times with IDD...", m, niter);
  else
    printf ("Multiplying %lu-bit number with %lu-bit number %lu times with IDD...", m, n, niter);
  fflush (stdout);
  
  reset_idd_stats();    
  t0 = cputime ();
  for (i = niter; i > 0; i--)
    z2 = x2*y2;    
  ti = cputime () - t0;
  printf ("done!\n");
  #ifdef MONITORING  
  reset_idd_stats();
  z2 = x2*y2;          
  printf("Before cleaning:\n");
  IDDContent::Ht.show_stats();  
  IDDContent::clean_Ht();
  printf("After cleaning:\n");
  IDDContent::Ht.show_stats();    
  show_idd_stats();
  #endif
  mul_per_sec_idd = 1000.0 * niter / ti;
  
  double vm2, res2;
  process_mem_usage(vm2,res2);

  printf ("RESULTS:\n");
  printf("x density: %f\n", (double)x2.a->density());
  printf("y density: %f\n", (double)y2.a->density());  
  printf("IDD: %f multiplication per second\n", mul_per_sec_idd);
  printf("Memory used: %dKo\n\n", (int)(res2-res1));
  #ifndef HIDE_RESULTS
  printf("x = "); mpz_out_str(stdout, 10, x);
  printf("\ny = "); mpz_out_str(stdout, 10, y);  
  printf("\nx*y = "); z2.a->print();
  printf("\n");    
  #endif
  mpz_mul(z,x,y);
  mpz_t comp;
  mpz_init(comp);
  z2.to_mpz(comp);
  if(mpz_cmp(comp,z)!=0) printf("Result is false!\n");
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
