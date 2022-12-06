#include<stdio.h>
#include<stdlib.h>

typedef struct 
{
  double r;
  double i;
} Complex;


Complex *
FFTmallocComplex1D(int num)
{
  Complex *ptr;
  /* malloc array  */
  if((ptr = (Complex *)malloc(sizeof(Complex) * num)) == NULL){
    fprintf(stderr, "Can't allocate memory.\n");
    exit(1);
  }
  return ptr;
} /* mallocDbl1D */


Complex ** 
FFTmallocComplex2D(int num1, int num2)
{
  Complex  *tmp1;
  Complex  **tmp2;
  int  i;

  tmp1 = (Complex *)malloc(sizeof(Complex) * num1 * num2);
  tmp2 = (Complex **)malloc(sizeof(Complex *) * num1);

  if((tmp1 == NULL)||(tmp2 == NULL)){
    fprintf(stderr, "Can't allocate memory.\n");
    exit(1);
  }

  for(i = 0 ; i < num1 ; i++)
    tmp2[i] = (Complex *)(tmp1 + i * num2);

  return(tmp2);
}

void
FFTfreeComplex1D(Complex *ptr)
{
  free(ptr);
}


void 
FFTfreeComplex2D(Complex **ptr)
{
  free(ptr[0]);
  free(ptr);
}



double *
FFTmallocDbl1D(int num)
{
  double *ptr;
  /* malloc array  */
  if((ptr = (double *)malloc(sizeof(double) * num)) == NULL){
    fprintf(stderr, "Can't allocate memory.\n");
    exit(1);
  }
  return ptr;
} /* mallocDbl1D */


double ** 
FFTmallocDbl2D(int num1, int num2)
{
  double  *tmp1;
  double  **tmp2;
  int  i;

  tmp1 = (double *)malloc(sizeof(double) * num1 * num2);
  tmp2 = (double **)malloc(sizeof(double *) * num1);

  if((tmp1 == NULL)||(tmp2 == NULL)){
    fprintf(stderr, "Can't allocate memory.\n");
    exit(1);
  }

  for(i = 0 ; i < num1 ; i++)
    tmp2[i] = (double *)(tmp1 + i * num2);

  return(tmp2);
}


void
FFTfreeDbl1D(double *ptr)
{
  free(ptr);
}


void 
FFTfreeDbl2D(double **ptr)
{
  free(ptr[0]);
  free(ptr);
}
