//
//  単母音認識プログラム
//  gcc -o vowel_recognition vowel_recognition.c wavheader.c -lm
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUF 256
#define N 32
#define PWR_THRES 500000
#define PWR_LAG 0.05
#define PATNUM 5

void readheader(int *samp, double *time, FILE *fp);

main(int argc,char *argv[]){
  FILE *fp, *fpv;  //ファイルポインタ
  // データ部
  char line[BUF];      // ファイルからのデータを読みこむ
  int samp;            // サンプリング周波数
  double time;         // 入力音声ファイルの時間長
  char vowel[] = {'a', 'i', 'u', 'e', 'o'};
  int i, j;
  char command[BUF];
  char patfname[BUF];
  double patdata[PATNUM][N];
  double testdata[N];
  double dist[PATNUM];         // ユークリッド距離
  double min_dist;
  int min_ind;
  char prefix[BUF];
  char infile[BUF];
  char envfile[BUF];
  double offset;
  double pwr, ppwr=0.0; 
  short s;

  if(argc != 3 && argc != 2){
    printf("usage: %s <inputfile prefix> [start time]\n", argv[0]);
    exit(1);
  }

  strcpy(prefix, argv[1]);
  sprintf(infile, "%s.wav", prefix);
  sprintf(envfile, "%s.env", prefix);

  // 音声開始位置検出
  if((fp=fopen(infile, "rb"))==NULL){
    printf("Can't open %s!\n", infile);
    exit(1);
  }
  readheader(&samp, &time, fp);
  //  samp = 16000;
  //  time = 0.6;
  i=0;
  while(i < samp * time){
    for(pwr = 0.0, j = 0; j < samp * PWR_LAG; i++, j++){
      fread(&s,sizeof(short),1,fp);
      pwr += s * s;
    }
    if(pwr - ppwr > PWR_THRES){
      printf("音声区間開始位置：　%.0lf ms\n", (double)i/samp*1000);
      offset = (double)i/samp;
      break;
    }
    ppwr = pwr;
  }
  if(argc==3)offset = atof(argv[2]);

  for(i = 0; i<PATNUM; i++){
    sprintf(patfname, "%c.env", vowel[i]);
    if((fpv=fopen(patfname, "r"))==NULL){
      printf("Can't open %s !\n", patfname);
      exit(1);
    }

    j=0;
    while(1){
      fgets(line, sizeof(line), fpv);
      if(feof(fpv))break;
      sscanf(line, "%lf", &patdata[i][j++]);
    }
    for(; j<N; j++){
      patdata[i][j] = 0.0;
    }
    fclose(fpv);
  }

  // ケプストラム分析
  sprintf(command, "./cepstrum.exe %s %lf", prefix, offset);
  system(command);

  // スペクトル包絡ファイルをオープン
  if((fp=fopen(envfile, "r"))==NULL){
    printf("Can't open %s !\n", envfile);
    exit(1);
  }
  
  j=0;
  while(1){
    fgets(line,sizeof(line),fp);
    if(feof(fp))break;
    sscanf(line, "%lf", &testdata[j++]);
  }
  for(; j<N; j++){
    testdata[j] = 0.0;
  }
  fclose(fp);

  for(i=0; i<PATNUM; i++){
    dist[i] = 0.0;
    for(j=0; patdata[i][j]>0 && testdata[j]>0; j++){
      dist[i] += (patdata[i][j] - testdata[j]) * (patdata[i][j] - testdata[j]);
    }
    dist[i] /= j;
    printf("[ %c ]との距離 = %lf\n", vowel[i], dist[i]);
  }

  for(min_dist=dist[0], min_ind=0, i=1; i<PATNUM; i++){
    if(min_dist>dist[i]){
      min_dist = dist[i];
      min_ind = i;
    }
  }

  // 結果を出力
  printf("\n認識結果：　[ %c ] (距離 = %lf)\n", vowel[min_ind], min_dist);

  fclose(fp);
 
  return 0;
}
