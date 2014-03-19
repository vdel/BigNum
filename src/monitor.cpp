#include "monitor.h"
#include <stdio.h>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>


int call_constructor, call_clean_Ht;
int call_T, call_int_T, call_COMP, call_C;
int call_OF, saved_call_OF;
int call_DEC, call_INC, call_int_INC;
int call_X, saved_call_X;
int call_RMSB, call_AMSB;
int call_TWICE, saved_call_TWICE;
int call_A, saved_call_A, call_int_A;
int call_P, saved_call_P, call_int_P;
int call_NEG;
int call_LSH, call_RSH;
int call_POW, call_Y, call_HALF;

void reset_idd_stats()
{
  call_constructor=0;
  call_clean_Ht=0;
  call_T=0;
  call_int_T=0;
  call_COMP=0;
  call_C=0;
  call_OF=0;
  saved_call_OF=0;
  call_DEC=0;
  call_INC=0;
  call_int_INC=0;
  call_X=0;
  saved_call_X=0;
  call_RMSB=0;
  call_AMSB=0;
  call_TWICE=0;
  saved_call_TWICE=0;
  call_A =0;
  call_int_A =0;  
  saved_call_A=0;
  call_P=0;
  call_int_P =0;  
  saved_call_P=0;
  call_NEG=0;
  call_LSH=0;
  call_RSH=0;
  call_POW =0;
  call_Y=0;
  call_HALF=0;
}

void show_idd_stats()
{
  printf("Stats:\n");
  printf("Constr: %d\n",call_constructor);
  printf("Clean Ht: %d\n",call_clean_Ht);
  printf("T: %d/%d\n",call_int_T,call_T);
  printf("C: %d\n",call_C);  
  printf("COMP: %d\n",call_COMP);
  printf("OF: %d/%d\n",call_OF,saved_call_OF);
  printf("DEC: %d\n",call_DEC);
  printf("INC: %d/%d\n",call_int_INC,call_INC);
  printf("X: %d/%d\n",call_X,saved_call_X);
  printf("RMSB: %d\n",call_RMSB);
  printf("AMSB: %d\n",call_AMSB);
  printf("TWICE: %d/%d\n",call_TWICE,saved_call_TWICE);
  printf("A: %d/%d/%d\n",call_int_A,call_A,saved_call_A);
  printf("P: %d/%d/%d\n",call_int_P,call_P,saved_call_P);
  printf("NEG: %d\n",call_NEG);
  printf("LSH: %d\n",call_LSH);
  printf("RSH: %d\n",call_RSH);
  printf("POW: %d\n",call_POW);
  printf("Y: %d\n",call_Y);
  printf("HALF: %d\n",call_HALF);
}

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned int vsize;
   int rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   int page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

