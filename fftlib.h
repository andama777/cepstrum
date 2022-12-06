/*==================================================================

  FFT library header file
    for fft1d.c, fft2d.c fftmem.c

  FFT library Lisence:      
    Copyright Takuya OOURA, 1996-1999 (fft4g.c fft4f2d.c)
    Modified by Hiroaki KAWASHIMA, 1999 (fft1d.c fft2d.c)
  
  Modification:
    Described in each source files.
==================================================================*/    

/* 複素数型を定義 */
typedef struct 
{
  double r;
  double i;
} Complex;


/*------------------------------------------------------------------
  １次元FFT, コサイン変換, サイン変換
------------------------------------------------------------------*/
void 
FFT1D(int, int, Complex *);

void 
DCT1D(int, int, double *);

void 
DST1D(int, int, double *);


/*------------------------------------------------------------------
  ２次元FFT, コサイン変換, サイン変換
------------------------------------------------------------------*/
void 
FFT2D(int, int, int, Complex **);

void 
DCT2D(int, int, int, double **);

void 
DST2D(int, int, int, double **);


/*------------------------------------------------------------------
  １、２次元メモリ配列確保

  FFTmallocComlex1D    複素数の１次元配列確保
  FFTmallocComlex2D    複素数の２次元配列確保
  FFTmallocDbl1D       実数の１次元配列確保（FFTmallocDbl1D 
  FFTmallocDbl2D       実数の２次元配列確保          ~~ l(エル) と 1(いち)

  FFTfreeComlex1D    複素数の１次元配列解放
  FFTfreeComlex2D    複素数の２次元配列解放
  FFTfreeDbl1D       実数の１次元配列解放（FFTfreeDbl1D 
  FFTfreeDbl2D       実数の２次元配列解放          ~~ l(エル) と 1(いち)

  使い方 : 例えば double a[n1][n2]; の代わりに
           double **a;
           a = FFTmallocDbl2D(n1, n2);
	   とする。

------------------------------------------------------------------*/

Complex *
FFTmallocComplex1D(int num);

Complex ** 
FFTmallocComplex2D(int num1, int num2);

void 
FFTfreeComplex1D(Complex *ptr);

void 
FFTfreeComplex2D(Complex **ptr);



double *
FFTmallocDbl1D(int num);

double ** 
FFTmallocDbl2D(int num1, int num2);

void 
FFTfreeDbl1D(double *ptr);

void 
FFTfreeDbl2D(double **ptr);



