/**
 * \file   mat-mult-pthrd.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for dense matrix multiplication using pthread.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-mult-pthrd.h"

#define __F4RT_DEBUG  0

#ifdef __F4RT_HAVE_PTHREAD_H

void *multPThreadImpose(void *p) {
  params *_p    = (params *)p;
  int tid       = _p->tid;
  uint32 start  = tid * _p->size;
  uint32 end    = start + _p->size;
  uint32 m      = _p->m;
  uint32 n      = _p->n;
  for (size_t i = start; i < end; ++i) {
    for(size_t j = 0; j < m; ++j) {
      mat sum = 0;
      for(size_t k = 0; k < n; ++k )
        sum += _p->a[k+i*n] * _p->b[k+j*n];
        //std::cout << j+i*m << "." << sum << std::endl;
      _p->c[j+i*m]  = sum;
    }
  }
  return 0;
}

void *multPThread(void *p) {
  params *_p    = (params *)p;
  int tid       = _p->tid;
  uint32 start  = tid * _p->size;
  uint32 end    = start + _p->size;
  uint32 m      = _p->m;
  uint32 n      = _p->n;
  for (size_t i = start; i < end; ++i) {
    for(size_t j = 0; j < m; ++j) {
      mat sum = 0;
      for(size_t k = 0; k < n; ++k )
        //std::cout << j+i*m << "." << sum << "LL" << std::endl;
        sum += _p->a[k+i*n] * _p->b[j+k*m];
      _p->c[j+i*m]  = sum;
    }
  }
  return 0;
}


// multiplies A*B^T and stores it in *this
void multPT(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
            uint32 blocksize, int impose) {
  uint32 l, m, n;
  if (impose == 1) {
    l = A.nRows();
    m = B.nRows();
    n = B.nCols();
  } else {
    l = A.nRows();
    m = B.nCols();
    n = B.nRows();
  }
  const mat *a_entries = A.entries.data();
  
  const mat *b_entries = new mat[n*m];
  b_entries = B.entries.data();

  mat *c_entries = new mat[l*m];

  std::cout << "Matrix Multiplication" << std::endl;
  timeval start, stop;
  clock_t cStart, cStop;
#if __F4RT_DEBUG
  std::cout << "A => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "B => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "C => " << C.nRows() << "-" << C.nCols() << "-" << C.nEntries() << std::endl;
#endif
  //const int padding = __F4RT_CPU_CACHE_LINE / sizeof(float);
  if (nthrds <= 0) {
    nthrds  = 1;
  }
  // holds thread information
  pthread_t threads[nthrds];
  params *thread_params = (params *) malloc(nthrds * sizeof(params));
  uint32 chunkSize  = l / nthrds;
  uint32 pad        = l % nthrds;

  // do matrix multiplication of submatrices of size in the order of 
  // blocksize
  gettimeofday(&start, NULL);
  cStart  = clock();
  if (impose == 1) {
    for (int i = 0; i < nthrds; ++i) {
      thread_params[i].a  = a_entries; 
      thread_params[i].b  = b_entries; 
      thread_params[i].c  = c_entries; 
      thread_params[i].tid  = i;
      // add 1 more chunk for the first pad threads
      if (i < pad)
        thread_params[i].size = chunkSize + 1;
      else
        thread_params[i].size = chunkSize;
      thread_params[i].m  = m;
      thread_params[i].n  = n;
      // real computation
      pthread_create(&threads[i], NULL, multPThreadImpose, (void *) &thread_params[i]);
    }

    // join threads back again
    for (int i = 0; i < nthrds; ++i)
      pthread_join(threads[i], NULL);

    free(thread_params);
  } else {
    for (int i = 0; i < nthrds; ++i) {
      thread_params[i].a  = a_entries; 
      thread_params[i].b  = b_entries; 
      thread_params[i].c  = c_entries; 
      thread_params[i].tid  = i;
      // add 1 more chunk for the first pad threads
      if (i < pad)
        thread_params[i].size = chunkSize + 1;
      else
        thread_params[i].size = chunkSize;
      thread_params[i].m  = m;
      thread_params[i].n  = n;
      // real computation
      pthread_create(&threads[i], NULL, multPThread, (void *) &thread_params[i]);
    }

    // join threads back again
    for (int i = 0; i < nthrds; ++i)
      pthread_join(threads[i], NULL);

    free(thread_params);
  }

  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           pThread 1D" << std::endl;
  std::cout << "Cache improved:   ";
  if (impose == 1)
    std::cout << "1" << std::endl;
  else
    std::cout << "0" << std::endl;
  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  float epsilon = 0.0000000001;
  double realtime = ((stop.tv_sec - start.tv_sec) * 1e6 + 
                    (stop.tv_usec - start.tv_usec)) / 1e6;
  double cputime  = (double)((cStop - cStart)) / CLOCKS_PER_SEC;
  char buffer[50];
  // get digits before decimal point of cputime (the longest number) and setw
  // with it: digits + 1 (point) + 4 (precision) 
  int digits = sprintf(buffer,"%.0f",cputime);
  double ratio = cputime/realtime;
  std::cout << "# Threads:        " << nthrds << std::endl;
  std::cout << "Block size:       " << blocksize << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "Real time:        " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << realtime << " sec" 
    << std::endl;
  std::cout << "CPU time:         " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << cputime
    << " sec" << std::endl;
  if (cputime > epsilon)
    std::cout << "CPU/real time:    " << std::setw(digits+1+4) 
      << std::setprecision(4) << std::fixed << ratio << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "GFLOPS/sec:       " << std::setw(digits+1+4) 
    << std::setprecision(4) << std::fixed << flops / (1000000000 * realtime) 
    << std:: endl;
  std::cout << "---------------------------------------------------" << std::endl;
}
#endif
