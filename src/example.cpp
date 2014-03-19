#include<stdio.h>
#include"iddlib.h"

uIDD ackermann(uIDD m, uIDD n);

int main(int, char**)
{
  uIDD fact = 1;
  for(int i=2;i<=1000;i++)
    fact = i*fact;
    
  printf("Factorielle 1000:\n"); idd_print(fact); printf("\n\n");
  
  uIDD fibo1 = 0;
  uIDD fibo2 = 1;
  for(int i=2;i<=10000;i++)
  {
    uIDD tmp(fibo2);
    fibo2 = fibo1+fibo2;
    fibo1 = tmp;
  }

  printf("Fibonacci 10000:\n"); idd_print(fibo2); printf("\n\n");
  
  printf("Ackermann(4,2):\n"); idd_print(ackermann(4,2)); printf("\n");  
  
  return 0;
}
  

uIDD ackermann(uIDD m, uIDD n)
{
  if(m == 0)
    return n+1;
  else if(m == 1)
    return n+2;
  else if(m == 2)
    return (2*n)+3;
  else if(m == 3)
    return (1<<(n+3))-3;
  else if(n == 0)
    return ackermann(m-1, 1);  
  else
    return ackermann(m-1, ackermann(m, n-1));
}
