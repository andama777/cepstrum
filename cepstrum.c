//
//  �P�v�X�g�������̓v���O����
//  gcc -o cepstrum cepstrum.c wavheader.c -lfftlib -lm -L.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fftlib.h"

#define BUF 256
#define LAG  512             // ���͑���
#define LAG2 32              // ��P�t�����V�p���i16000Hz/400Hz��2�̋Џ�ɋ߂��l�j
#define PI 3.14159
#define MAX_F0 400           // ��{���g���̍ő�l�i���o����臒l�j

void readheader(int *samp, double *time, FILE *fp);

main(int argc, char *argv[]){
  FILE *fp, *fp1, *fp2,*fp3;//�t�@�C���|�C���^
  char infile[BUF];     //���̓t�@�C����
  char outfile1[BUF];   //�o�̓t�@�C��1
  char outfile2[BUF];   //�o�̓t�@�C��2
  char outfile3[BUF];   //�o�̓t�@�C��3
  // �w�b�_��񂩂� 
  int samp;            //�T���v�����O���g��   
  double time;          //����
  //  �v���O�������f�[�^
  Complex *c1,*c2,*c3,*c4,*c5;  // �����f�[�^(���֐��A�t�[���G�ϊ��p)
  short *sfp;               // �����f�[�^�܂邲��         
  int i,j,k;                // ���[�v�Q�Ɨp�֐�
  double max = -10000.0;           // ���P�t�����V
  int max_label;            // ��Ԓ�P�t�����V�̂Ƃ��̃��x��
  double offset = 0.0;

  if(argc != 2 && argc != 3){
    printf("usage: %s <input wavfile> [offset]\n", argv[0]);
    exit(1);
  }

  sprintf(infile, "%s.wav", argv[1]);
  sprintf(outfile1, "%s.env", argv[1]);
  sprintf(outfile2, "%s.f0", argv[1]);
  sprintf(outfile3, "%s.cep", argv[1]);
  if(argc==3)offset = atof(argv[2]);

  // �����t�@�C�����I�[�v��
  
  if((fp=fopen(infile, "rb"))==NULL){
    printf("Can't open %s !\n", infile);
    exit(1);
  }
  
  // �X�y�N�g����������o���t�@�C�����J��
  if((fp1=fopen(outfile1, "w"))==NULL){
    printf("Can't open %s !\n", outfile1);
    exit(1);
  }
  
  // ��{���g���������o���t�@�C�����J��
  if((fp2=fopen(outfile2, "w"))==NULL){
    printf("Can't open %s !\n", outfile2);
    exit(1);
  }

  // �P�t�����V�f�[�^�������o���t�@�C�����J��
  if((fp3=fopen(outfile3, "w"))==NULL){
    printf("Can't open %s !\n", outfile3);
    exit(1);
  }

  readheader(&samp, &time, fp);
  //  samp = 16000;
  //  time = 0.6;

  // �����f�[�^���ׂĂ��i�[���镪�������I�ɔz�������
  sfp = (short*)malloc(sizeof(short)*(samp*time));

  // �t�[���G�ϊ�����f�[�^�̕��������I�ɔz�������
  c1 = FFTmallocComplex1D(LAG);    // ������
  c2 = FFTmallocComplex1D(LAG);    // DFT
  c3 = FFTmallocComplex1D(LAG);    // �ΐ��ϊ��{IDFT
  c4 = FFTmallocComplex1D(LAG);    // �X�y�N�g���DFT�v�Z�p
  c5 = FFTmallocComplex1D(LAG);    // �w���ϊ��̌���(�X�y�N�g���)

   // �����f�[�^��z��ɂ����
  for(i = 0; i < samp * time; i++){
    fread(sfp+i,sizeof(short),1,fp); 
  }
  
  // �������͋�Ԃ̎n�_�Ƀ|�C���^���ړ�
  sfp += (int)(offset * samp);  
 
  // ��������Z���ԃX�y�N�g�����́B���P�t�����V�����߂�
  // �t���[���̒������͂��ׂ�1�̏����ōs�Ȃ���̂ł��̕�(LAG)
  for(j = 0; j < LAG; j++){
    // ��������͑��֐�(�n�~���O��)
    c1[j].r = *(sfp+j) * (0.54 - 0.46 * cos(2*PI*j/(LAG-1)));
    c2[j].r = c1[j].r;
    c2[j].i = c1[j].i = 0.0;
  }

  // ��������DFT
  FFT1D(LAG, -1, c2);
 
  // ��������ΐ��ϊ� (��Βl�̑ΐ�)��Βl��r,i�̂Q��̘a�́���Ƃ�
  for(j = 0; j < LAG; j++){
    c3[j].r = log(sqrt(c2[j].r * c2[j].r + c2[j].i * c2[j].i));
    c3[j].i = 0;
  }     
 
  // ��������IDFT
  FFT1D(LAG, 1, c3);
  
  for(j = 0; j < LAG; j++){
    c4[j].i = 0.0;                        //
    c4[j].r = c3[j].r /= (LAG);           // IDFT�Ȃ̂� 1/N ��������
    fprintf(fp3,"%lf\n",c3[j].r);
  }   
    
  //  �s�[�N���o
  max_label = samp/MAX_F0;
  for(j = samp/MAX_F0; j < LAG/2; j++){
    if(max < c3[j].r){
      max = c3[j].r;
      max_label = j;
    }
  }  
  // ��{���g�����o��
  fprintf(fp2,"F0 = %d Hz(= %d Hz/ %d sample point)\n",
	  samp/max_label, samp, max_label);
  printf("\nF0 = %d Hz(= %d Hz/ %d sample point)\n",
	 samp/max_label, samp, max_label);

  // �s�[�N���o�����܂� 

  // ��������X�y�N�g������o
  // ��������͑��֐�(�n�~���O��)(���t�^�����O)  
  for(j = 0; j < LAG2; j++){
    c4[j].r = c4[j].r * (0.54 - 0.46 * cos(2*PI*j/(LAG-1)));
  }

  // �t�[���G�ϊ����s�Ȃ�
  FFT1D(LAG2, -1, c4);

  // �w���ϊ����������ʂ��X�y�N�g����̌��ʂƂȂ�
  for(j = 0; j < LAG2/2; j++){
    c5[j].r = exp(c4[j].r);
    //    printf("spectol[%d]==>%lf\n",j,c5[j].r);
    fprintf(fp1,"%lf\n",c5[j].r);
  }  
  // �X�y�N�g������o�����܂�

  FFTfreeComplex1D(c1);
  FFTfreeComplex1D(c2);
  FFTfreeComplex1D(c3);
  FFTfreeComplex1D(c4);
  FFTfreeComplex1D(c5);
 
  fclose(fp);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

  return 0;
}
