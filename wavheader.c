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
  
  ////////// �w�b�_ ////////////
  fwrite(s_riff, sizeof(char), 4, fp); // "RIFF"

  d = samp*time*2+36; // �ȍ~�̃f�[�^�T�C�Y�i�������f�[�^�T�C�Y�{�w�b�_�̎c��36Bytes�j
  fwrite(&d, sizeof(int), 1, fp);
  
  ////////// �����f�[�^��� ////////////
  fwrite(s_wave, sizeof(char), 8, fp); // "WAVEfmt "

  d = 16; // �����f�[�^���̃o�C�g��
  fwrite(&d, sizeof(int), 1, fp);

  s = 1;  // �f�[�^�`���iPCM��1�j
  fwrite(&s, sizeof(short), 1, fp);

  s = 1;  // �`�����l�����i���m������1�A�X�e���I��2�j
  fwrite(&s, sizeof(short), 1, fp);

  d = samp;  // �T���v�����O���g���iHz�j
  fwrite(&d, sizeof(int), 1, fp);

  d = samp*2;  // 1�b������̉����f�[�^�̃o�C�g���i16bit�Ȃ̂Ł~2�j
  fwrite(&d, sizeof(int), 1, fp);

  s = 2;  // 1�T���v��������̃o�C�g���i16bit�Ȃ̂�2�j
  fwrite(&s, sizeof(short), 1, fp);

  s = 16;  // �ʎq���r�b�g��
  fwrite(&s, sizeof(short), 1, fp);

  ////////// �����f�[�^ ////////////
  fwrite(s_data, sizeof(char), 4, fp); // "data"

  d = samp*time*2;  // �����f�[�^�̃o�C�g��
  fwrite(&d, sizeof(int), 1, fp);
}

void readheader(int *samp, double *time, FILE *fp){
  char buf[BUF];
  short s;
  int d, channel;
  
  ////////// �w�b�_ ////////////
  fread(buf, sizeof(char), 4, fp); // "RIFF"
  buf[4] = '\0';
  //  printf("�����w�b�_����\n��񎯕ʎq: [%s]\n", buf);

  //  d = samp*time*2+36; // �ȍ~�̃f�[�^�T�C�Y�i�������f�[�^�T�C�Y�{�w�b�_�̎c��36Bytes�j
  fread(&d, sizeof(int), 1, fp);
  //printf("�ȍ~�̃f�[�^�T�C�Y: [%d]�o�C�g�i�������f�[�^�T�C�Y�{�w�b�_�̎c��36Bytes�j\n", d);

  ////////// �����f�[�^��� ////////////
  fread(buf, sizeof(char), 8, fp); // "WAVEfmt "
  buf[8] = '\0';
  //printf("�����f�[�^�t�H�[�}�b�g: [%s]\n", buf);

  //d = 16; // �����f�[�^���̃o�C�g��
  fread(&d, sizeof(int), 1, fp);
  //printf("�����f�[�^���̃o�C�g��: [%d]�o�C�g\n", d);
  
  //s = 1;  // �f�[�^�`���iPCM��1�j
  fread(&s, sizeof(short), 1, fp);
  //printf("�f�[�^�`��: [%d]�i��PCM�j\n", s);

  //  s = 1;  // �`�����l�����i���m������1�A�X�e���I��2�j
  fread(&s, sizeof(short), 1, fp);
  //printf("�`�����l����: [%d]�i���m������1�A�X�e���I��2�j\n", s);
  channel = (int)s;

  //d = samp;  // �T���v�����O���g���iHz�j
  fread(&d, sizeof(int), 1, fp);
  //printf("�T���v�����O���g��: [%d] �iHz�j\n", d);
  *samp = d;

  //d = samp*2;  // 1�b������̉����f�[�^�̃o�C�g���i16bit�Ȃ̂Ł~2�j
  fread(&d, sizeof(int), 1, fp);
  //printf("1�b������̉����f�[�^�̃o�C�g��: [%d]�o�C�g\n", d);

  //  s = 2;  // 1�T���v��������̃o�C�g���i16bit�Ȃ̂�2�j
  fread(&s, sizeof(short), 1, fp);
  //printf("1�T���v��������̃o�C�g��: [%d]�o�C�g�i��16bit�j\n", s);

  // s = 16;  // �ʎq���r�b�g��
  fread(&s, sizeof(short), 1, fp);
  //printf("�ʎq���r�b�g��: [%d] bit\n", s);

  ////////// �����f�[�^ ////////////
  fread(buf, sizeof(char), 4, fp); // "data"
  buf[4] = '\0';
  //printf("�f�[�^���ʎq: [%s]\n", buf);
  
  //d = samp*time*2;  // �����f�[�^�̃o�C�g��
  fread(&d, sizeof(int), 1, fp);
  //printf("�����f�[�^�̃o�C�g��: [%d]�o�C�g\n", d);
  *time = (double)d / *samp / channel / 2;
}

void readheader2(int *samp, int *size, FILE *fp){
  char buf[BUF];
  short s;
  int d, channel;
  
  ////////// �w�b�_ ////////////
  fread(buf, sizeof(char), 4, fp); // "RIFF"
  buf[4] = '\0';
  printf("�����w�b�_����\n��񎯕ʎq: [%s]\n", buf);

  //  d = samp*time*2+36; // �ȍ~�̃f�[�^�T�C�Y�i�������f�[�^�T�C�Y�{�w�b�_�̎c��36Bytes�j
  fread(&d, sizeof(int), 1, fp);
  printf("�ȍ~�̃f�[�^�T�C�Y: [%d]�o�C�g�i�������f�[�^�T�C�Y�{�w�b�_�̎c��36Bytes�j\n", d);

  ////////// �����f�[�^��� ////////////
  fread(buf, sizeof(char), 8, fp); // "WAVEfmt "
  buf[8] = '\0';
  printf("�����f�[�^�t�H�[�}�b�g: [%s]\n", buf);

  //d = 16; // �����f�[�^���̃o�C�g��
  fread(&d, sizeof(int), 1, fp);
  printf("�����f�[�^���̃o�C�g��: [%d]�o�C�g\n", d);
  
  //s = 1;  // �f�[�^�`���iPCM��1�j
  fread(&s, sizeof(short), 1, fp);
  printf("�f�[�^�`��: [%d]�i��PCM�j\n", s);

  //  s = 1;  // �`�����l�����i���m������1�A�X�e���I��2�j
  fread(&s, sizeof(short), 1, fp);
  printf("�`�����l����: [%d]�i���m������1�A�X�e���I��2�j\n", s);
  channel = (int)s;

  //d = samp;  // �T���v�����O���g���iHz�j
  fread(&d, sizeof(int), 1, fp);
  printf("�T���v�����O���g��: [%d] �iHz�j\n", d);
  *samp = d;

  //d = samp*2;  // 1�b������̉����f�[�^�̃o�C�g���i16bit�Ȃ̂Ł~2�j
  fread(&d, sizeof(int), 1, fp);
  printf("1�b������̉����f�[�^�̃o�C�g��: [%d]�o�C�g\n", d);

  //  s = 2;  // 1�T���v��������̃o�C�g���i16bit�Ȃ̂�2�j
  fread(&s, sizeof(short), 1, fp);
  printf("1�T���v��������̃o�C�g��: [%d]�o�C�g�i��16bit�j\n", s);

  // s = 16;  // �ʎq���r�b�g��
  fread(&s, sizeof(short), 1, fp);
  printf("�ʎq���r�b�g��: [%d] bit\n", s);

  ////////// �����f�[�^ ////////////
  fread(buf, sizeof(char), 4, fp); // "data"
  buf[4] = '\0';
  printf("�f�[�^���ʎq: [%s]\n", buf);
  
  //d = samp*time*2;  // �����f�[�^�̃o�C�g��
  fread(&d, sizeof(int), 1, fp);
  printf("�����f�[�^�̃o�C�g��: [%d]�o�C�g\n", d);
  *size = d / channel / 2;
}
