#include <stdio.h>
#include "cs.h"

// make
// gcc -Wall -O test_csparse.c -o test_csparse -L. -lcsparse
// test_csparse < t1

int main ()
{
  cs *M, *A;
  csd *BM;

  M = cs_load(stdin);

  cs_print(M,1);

  A = cs_compress(M);

  BM = cs_dmperm(A,1);
  //BM = cs_scc(A);

  if (!BM) { fprintf(stderr,"Error\n"); exit(1); }

  printf("nb = %d\n",(int)BM->nb);

  int i;
  for (i = 1; i <= BM->nb; ++i)
    {
      printf("%d ",(int)BM->r[i]);
    }

  printf("\n");

  return 0;
}
