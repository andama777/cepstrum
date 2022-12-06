//
//  ケプストラム分析プログラム
//  gcc -o cepstrum cepstrum.c wavheader.c -lfftlib -lm -L.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fftlib.h"

#define BUF 256
#define LAG  512             // 分析窓長
#define LAG2 32              // 低ケフランシ用窓（16000Hz/400Hzで2の巾乗に近い値）
#define PI 3.14159
#define MAX_F0 400           // 基本周波数の最大値（検出時の閾値）

void readheader(int *samp, double *time, FILE *fp);

main(int argc, char *argv[]){
  FILE *fp, *fp1, *fp2,*fp3;//ファイルポインタ
  char infile[BUF];     //入力ファイル名
  char outfile1[BUF];   //出力ファイル1
  char outfile2[BUF];   //出力ファイル2
  char outfile3[BUF];   //出力ファイル3
  // ヘッダ情報から 
  int samp;            //サンプリング周波数   
  double time;          //時間
  //  プログラム内データ
  Complex *c1,*c2,*c3,*c4,*c5;  // 音声データ(窓関数、フーリエ変換用)
  short *sfp;               // 音声データまるごと         
  int i,j,k;                // ループ参照用関数
  double max = -10000.0;           // 高ケフランシ
  int max_label;            // 一番低ケフランシのとこのラベル
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

  // 音声ファイルをオープン
  
  if((fp=fopen(infile, "rb"))==NULL){
    printf("Can't open %s !\n", infile);
    exit(1);
  }
  
  // スペクトル包絡を書き出すファイルを開く
  if((fp1=fopen(outfile1, "w"))==NULL){
    printf("Can't open %s !\n", outfile1);
    exit(1);
  }
  
  // 基本周波数を書き出すファイルを開く
  if((fp2=fopen(outfile2, "w"))==NULL){
    printf("Can't open %s !\n", outfile2);
    exit(1);
  }

  // ケフランシデータを書き出すファイルを開く
  if((fp3=fopen(outfile3, "w"))==NULL){
    printf("Can't open %s !\n", outfile3);
    exit(1);
  }

  readheader(&samp, &time, fp);
  //  samp = 16000;
  //  time = 0.6;

  // 音声データすべてを格納する分だけ動的に配列をつくる
  sfp = (short*)malloc(sizeof(short)*(samp*time));

  // フーリエ変換するデータの分だけ動的に配列をつくる
  c1 = FFTmallocComplex1D(LAG);    // 窓がけ
  c2 = FFTmallocComplex1D(LAG);    // DFT
  c3 = FFTmallocComplex1D(LAG);    // 対数変換＋IDFT
  c4 = FFTmallocComplex1D(LAG);    // スペクトル包絡DFT計算用
  c5 = FFTmallocComplex1D(LAG);    // 指数変換の結果(スペクトル包絡)

   // 音声データを配列にいれる
  for(i = 0; i < samp * time; i++){
    fread(sfp+i,sizeof(short),1,fp); 
  }
  
  // 音声分析区間の始点にポインタを移動
  sfp += (int)(offset * samp);  
 
  // ここから短時間スペクトル分析。高ケフランシを求める
  // フレームの長さ分はすべて1つの処理で行なえるのでその分(LAG)
  for(j = 0; j < LAG; j++){
    // ここからは窓関数(ハミング窓)
    c1[j].r = *(sfp+j) * (0.54 - 0.46 * cos(2*PI*j/(LAG-1)));
    c2[j].r = c1[j].r;
    c2[j].i = c1[j].i = 0.0;
  }

  // ここからDFT
  FFT1D(LAG, -1, c2);
 
  // ここから対数変換 (絶対値の対数)絶対値はr,iの２乗の和の√をとる
  for(j = 0; j < LAG; j++){
    c3[j].r = log(sqrt(c2[j].r * c2[j].r + c2[j].i * c2[j].i));
    c3[j].i = 0;
  }     
 
  // ここからIDFT
  FFT1D(LAG, 1, c3);
  
  for(j = 0; j < LAG; j++){
    c4[j].i = 0.0;                        //
    c4[j].r = c3[j].r /= (LAG);           // IDFTなので 1/N をかける
    fprintf(fp3,"%lf\n",c3[j].r);
  }   
    
  //  ピーク検出
  max_label = samp/MAX_F0;
  for(j = samp/MAX_F0; j < LAG/2; j++){
    if(max < c3[j].r){
      max = c3[j].r;
      max_label = j;
    }
  }  
  // 基本周波数を出力
  fprintf(fp2,"F0 = %d Hz(= %d Hz/ %d sample point)\n",
	  samp/max_label, samp, max_label);
  printf("\nF0 = %d Hz(= %d Hz/ %d sample point)\n",
	 samp/max_label, samp, max_label);

  // ピーク抽出おしまい 

  // ここからスペクトル包絡抽出
  // ここからは窓関数(ハミング窓)(リフタリング)  
  for(j = 0; j < LAG2; j++){
    c4[j].r = c4[j].r * (0.54 - 0.46 * cos(2*PI*j/(LAG-1)));
  }

  // フーリエ変換を行なう
  FFT1D(LAG2, -1, c4);

  // 指数変換をした結果がスペクトル包絡の結果となる
  for(j = 0; j < LAG2/2; j++){
    c5[j].r = exp(c4[j].r);
    //    printf("spectol[%d]==>%lf\n",j,c5[j].r);
    fprintf(fp1,"%lf\n",c5[j].r);
  }  
  // スペクトル包絡抽出おしまい

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
