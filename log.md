# until 00:00

```
make
gcc -o cepstrum cepstrum.c wavheader.c -lfftlib -lm -L.
gcc -o vowel_recognition vowel_recognition.c wavheader.c -lm
```