#ifndef F4RT_MAT_MULT_PTHRD_H
#define F4RT_MAT_MULT_PTHRD_H

#include <matrix.h>

struct params {
  const float *a;
  const float *b;
  float *c;
  int size;
  uint32 m;
  uint32 n;
  int tid;
};

void *multPThreadImpose(void *p);

void *multPThread(void *p);

// multiplies A*B^T and stores it in C
void multPT(Matrix& C, const Matrix& A, const Matrix& B, int nthrds,
            int blocksize, int impose);

#endif
