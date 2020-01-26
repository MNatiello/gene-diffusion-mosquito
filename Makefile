INCLUDE=../dSFMT-src-2.2.3
oxiV5:
	gcc -O3 -o Oxi-FellerV5 -lm -IINCLUDE -LINCLUDE -DDSFMT_MEXP=19937  Oxi-FellerV5.c ${INCLUDE}/dSFMT.o
