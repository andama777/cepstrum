#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUF 256

void makeheader(int samp, double time, FILE *fp){
  char s_riff[] = "RIFF";
  char s_wave[] = "WAVEfmt ";
  char s_data[] = "data";
  short s;
  int d;
  
  ////////// ヘッダ ////////////
  fwrite(s_riff, sizeof(char), 4, fp); // "RIFF"

  d = samp*time*2+36; // 以降のデータサイズ（＝音声データサイズ＋ヘッダの残り36Bytes）
  fwrite(&d, sizeof(int), 1, fp);
  
  ////////// 音声データ情報 ////////////
  fwrite(s_wave, sizeof(char), 8, fp); // "WAVEfmt "

  d = 16; // 音声データ情報のバイト数
  fwrite(&d, sizeof(int), 1, fp);

  s = 1;  // データ形式（PCM＝1）
  fwrite(&s, sizeof(short), 1, fp);

  s = 1;  // チャンネル数（モノラル＝1、ステレオ＝2）
  fwrite(&s, sizeof(short), 1, fp);

  d = samp;  // サンプリング周波数（Hz）
  fwrite(&d, sizeof(int), 1, fp);

  d = samp*2;  // 1秒あたりの音声データのバイト数（16bitなので×2）
  fwrite(&d, sizeof(int), 1, fp);

  s = 2;  // 1サンプルあたりのバイト数（16bitなので2）
  fwrite(&s, sizeof(short), 1, fp);

  s = 16;  // 量子化ビット数
  fwrite(&s, sizeof(short), 1, fp);

  ////////// 音声データ ////////////
  fwrite(s_data, sizeof(char), 4, fp); // "data"

  d = samp*time*2;  // 音声データのバイト数
  fwrite(&d, sizeof(int), 1, fp);
}

void readheader(int *samp, double *time, FILE *fp){
  char buf[BUF];
  short s;
  int d, channel;
  
  ////////// ヘッダ ////////////
  fread(buf, sizeof(char), 4, fp); // "RIFF"
  buf[4] = '\0';
  //  printf("＜＜ヘッダ＞＞\n情報識別子: [%s]\n", buf);

  //  d = samp*time*2+36; // 以降のデータサイズ（＝音声データサイズ＋ヘッダの残り36Bytes）
  fread(&d, sizeof(int), 1, fp);
  //printf("以降のデータサイズ: [%d]バイト（＝音声データサイズ＋ヘッダの残り36Bytes）\n", d);

  ////////// 音声データ情報 ////////////
  fread(buf, sizeof(char), 8, fp); // "WAVEfmt "
  buf[8] = '\0';
  //printf("音声データフォーマット: [%s]\n", buf);

  //d = 16; // 音声データ情報のバイト数
  fread(&d, sizeof(int), 1, fp);
  //printf("音声データ情報のバイト数: [%d]バイト\n", d);
  
  //s = 1;  // データ形式（PCM＝1）
  fread(&s, sizeof(short), 1, fp);
  //printf("データ形式: [%d]（＝PCM）\n", s);

  //  s = 1;  // チャンネル数（モノラル＝1、ステレオ＝2）
  fread(&s, sizeof(short), 1, fp);
  //printf("チャンネル数: [%d]（モノラル＝1、ステレオ＝2）\n", s);
  channel = (int)s;

  //d = samp;  // サンプリング周波数（Hz）
  fread(&d, sizeof(int), 1, fp);
  //printf("サンプリング周波数: [%d] （Hz）\n", d);
  *samp = d;

  //d = samp*2;  // 1秒あたりの音声データのバイト数（16bitなので×2）
  fread(&d, sizeof(int), 1, fp);
  //printf("1秒あたりの音声データのバイト数: [%d]バイト\n", d);

  //  s = 2;  // 1サンプルあたりのバイト数（16bitなので2）
  fread(&s, sizeof(short), 1, fp);
  //printf("1サンプルあたりのバイト数: [%d]バイト（＝16bit）\n", s);

  // s = 16;  // 量子化ビット数
  fread(&s, sizeof(short), 1, fp);
  //printf("量子化ビット数: [%d] bit\n", s);

  ////////// 音声データ ////////////
  fread(buf, sizeof(char), 4, fp); // "data"
  buf[4] = '\0';
  //printf("データ識別子: [%s]\n", buf);
  
  //d = samp*time*2;  // 音声データのバイト数
  fread(&d, sizeof(int), 1, fp);
  //printf("音声データのバイト数: [%d]バイト\n", d);
  *time = (double)d / *samp / channel / 2;
}

void readheader2(int *samp, int *size, FILE *fp){
  char buf[BUF];
  short s;
  int d, channel;
  
  ////////// ヘッダ ////////////
  fread(buf, sizeof(char), 4, fp); // "RIFF"
  buf[4] = '\0';
  printf("＜＜ヘッダ＞＞\n情報識別子: [%s]\n", buf);

  //  d = samp*time*2+36; // 以降のデータサイズ（＝音声データサイズ＋ヘッダの残り36Bytes）
  fread(&d, sizeof(int), 1, fp);
  printf("以降のデータサイズ: [%d]バイト（＝音声データサイズ＋ヘッダの残り36Bytes）\n", d);

  ////////// 音声データ情報 ////////////
  fread(buf, sizeof(char), 8, fp); // "WAVEfmt "
  buf[8] = '\0';
  printf("音声データフォーマット: [%s]\n", buf);

  //d = 16; // 音声データ情報のバイト数
  fread(&d, sizeof(int), 1, fp);
  printf("音声データ情報のバイト数: [%d]バイト\n", d);
  
  //s = 1;  // データ形式（PCM＝1）
  fread(&s, sizeof(short), 1, fp);
  printf("データ形式: [%d]（＝PCM）\n", s);

  //  s = 1;  // チャンネル数（モノラル＝1、ステレオ＝2）
  fread(&s, sizeof(short), 1, fp);
  printf("チャンネル数: [%d]（モノラル＝1、ステレオ＝2）\n", s);
  channel = (int)s;

  //d = samp;  // サンプリング周波数（Hz）
  fread(&d, sizeof(int), 1, fp);
  printf("サンプリング周波数: [%d] （Hz）\n", d);
  *samp = d;

  //d = samp*2;  // 1秒あたりの音声データのバイト数（16bitなので×2）
  fread(&d, sizeof(int), 1, fp);
  printf("1秒あたりの音声データのバイト数: [%d]バイト\n", d);

  //  s = 2;  // 1サンプルあたりのバイト数（16bitなので2）
  fread(&s, sizeof(short), 1, fp);
  printf("1サンプルあたりのバイト数: [%d]バイト（＝16bit）\n", s);

  // s = 16;  // 量子化ビット数
  fread(&s, sizeof(short), 1, fp);
  printf("量子化ビット数: [%d] bit\n", s);

  ////////// 音声データ ////////////
  fread(buf, sizeof(char), 4, fp); // "data"
  buf[4] = '\0';
  printf("データ識別子: [%s]\n", buf);
  
  //d = samp*time*2;  // 音声データのバイト数
  fread(&d, sizeof(int), 1, fp);
  printf("音声データのバイト数: [%d]バイト\n", d);
  *size = d / channel / 2;
}
