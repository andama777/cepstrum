/*==================================================================

  FFT library header file
    for fft1d.c, fft2d.c fftmem.c

  FFT library Lisence:      
    Copyright Takuya OOURA, 1996-1999 (fft4g.c fft4f2d.c)
    Modified by Hiroaki KAWASHIMA, 1999 (fft1d.c fft2d.c)
  
  Modification:
    Described in each source files.
==================================================================*/    

/* $BJ#AG?t7?$rDj5A(B */
typedef struct 
{
  double r;
  double i;
} Complex;


/*------------------------------------------------------------------
  $B#1<!85(BFFT, $B%3%5%$%sJQ49(B, $B%5%$%sJQ49(B
------------------------------------------------------------------*/
void 
FFT1D(int, int, Complex *);

void 
DCT1D(int, int, double *);

void 
DST1D(int, int, double *);


/*------------------------------------------------------------------
  $B#2<!85(BFFT, $B%3%5%$%sJQ49(B, $B%5%$%sJQ49(B
------------------------------------------------------------------*/
void 
FFT2D(int, int, int, Complex **);

void 
DCT2D(int, int, int, double **);

void 
DST2D(int, int, int, double **);


/*------------------------------------------------------------------
  $B#1!"#2<!85%a%b%jG[Ns3NJ](B

  FFTmallocComlex1D    $BJ#AG?t$N#1<!85G[Ns3NJ](B
  FFTmallocComlex2D    $BJ#AG?t$N#2<!85G[Ns3NJ](B
  FFTmallocDbl1D       $B<B?t$N#1<!85G[Ns3NJ]!J(BFFTmallocDbl1D 
  FFTmallocDbl2D       $B<B?t$N#2<!85G[Ns3NJ](B          ~~ l($B%(%k(B) $B$H(B 1($B$$$A(B)

  FFTfreeComlex1D    $BJ#AG?t$N#1<!85G[Ns2rJ|(B
  FFTfreeComlex2D    $BJ#AG?t$N#2<!85G[Ns2rJ|(B
  FFTfreeDbl1D       $B<B?t$N#1<!85G[Ns2rJ|!J(BFFTfreeDbl1D 
  FFTfreeDbl2D       $B<B?t$N#2<!85G[Ns2rJ|(B          ~~ l($B%(%k(B) $B$H(B 1($B$$$A(B)

  $B;H$$J}(B : $BNc$($P(B double a[n1][n2]; $B$NBe$o$j$K(B
           double **a;
           a = FFTmallocDbl2D(n1, n2);
	   $B$H$9$k!#(B

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



