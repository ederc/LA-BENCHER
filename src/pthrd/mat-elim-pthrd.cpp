/**
 * \file   mat-elim-pthrd.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for Gaussian Elimination using pThread.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "mat-elim-pthrd.h"

#define F4RT_DBG  0

#ifdef __F4RT_HAVE_PTHREAD_H


// counter for number of active threads
pthread_mutex_t mutex1  = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond1    = PTHREAD_COND_INITIALIZER;
int counter             = 0;
int jobsRunning         = 0;
int maxThreads          = 0;

void *elimPTHRD(void *p) {
  paramsElim *_p  = (paramsElim *)p;
  uint32 start    = _p->start + _p->index;
  uint32 end      = start + _p->size;
  uint32 n        = _p->n;
  uint32 i        = _p->index;
  uint64 prime    = _p->prime;
  mat inv         = _p->inv;
  mat mult;
  for (uint32 j = start+1; j < end+1; ++j) {
    mult  = (_p->a[i+j*n] * inv) % prime;
    for (uint32 k = i+1; k < n; ++k) {
      _p->a[k+j*n]  +=  _p->a[k+i*n] * mult;
      _p->a[k+j*n]  %=  prime;
    }
  }
  return 0;
}

void elimNaivePTHRDModP1d(Matrix& A, int numberThreads, uint32 blocksize, uint64 prime) {
  uint32 nthrds   = numberThreads;
  uint32 m        = A.nRows();
  uint32 n        = A.nCols();
  mat *a_entries  = A.entries.data();
  // if m > n then only n eliminations are possible
  uint32 boundary  = (m > n) ? n : m;
  mat inv;
  //const int padding = __F4RT_CPU_CACHE_LINE / sizeof(float);
  if (nthrds <= 0) {
    nthrds  = 1;
  }
  // holds thread information
  pthread_t threads[nthrds];
  paramsElim *thread_params = (paramsElim *) malloc(nthrds * sizeof(paramsElim));
  uint32 chunkSize;
  uint32 pad;
  uint32 ctr;

  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Naive Gaussian Elimination without pivoting" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();
  for (uint32 i = 0; i < boundary; ++i) {
    A(i,i) %= prime;
#if F4RT_DBG
    std::cout << "!! A(" << i << "," << i << ") " << A(i,i) << std::endl;
    std::cout << "A(" << i << "," << i << ") " << A(i,i) % prime << std::endl;
#endif
    inv  = negInverseModP(A(i,i), prime);
#if F4RT_DBG
    std::cout << "inv  " << inv << std::endl;
#endif
    chunkSize = (m - i - 1) / nthrds;
    pad       = (m - i - 1) % nthrds;
    ctr = 0;
    for (uint32 l = 0; l < nthrds; ++l) {
      thread_params[l].a      = a_entries;
      thread_params[l].prime  = prime;
      thread_params[l].index  = i;
      if (l < pad) {
        thread_params[l].size   = chunkSize + 1;
        thread_params[l].start  = ctr;
        ctr +=  chunkSize + 1;
      } else {
        thread_params[l].size = chunkSize;
        thread_params[l].start  = ctr;
        ctr +=  chunkSize;
      }
      thread_params[l].n    = n;
      thread_params[l].inv  = inv;
      // real computation
      pthread_create(&threads[l], NULL, elimPTHRD, (void *) &thread_params[l]);
    }

    // join threads back again
    for (uint32 l = 0; l < nthrds; ++l)
      pthread_join(threads[l], NULL);
  }
  free(thread_params);
  //cleanUpModP(A, prime);
  //A.print();
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           pThread 1D" << std::endl;
  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = countGEPFlops(m, n, prime);
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

void elimNaivePTHRDModP1dPivot(Matrix& A, int numberThreads, uint32 blocksize, uint64 prime) {
  uint32 nthrds   = numberThreads;
  uint32 l;
  uint32 m        = A.nRows();
  uint32 n        = A.nCols();
  mat *a_entries  = A.entries.data();
  // if m > n then only n eliminations are possible
  uint32 boundary  = (m > n) ? n : m;
  mat inv;
  //const int padding = __F4RT_CPU_CACHE_LINE / sizeof(float);
  if (nthrds <= 0) {
    nthrds  = 1;
  }
  // holds thread information
  pthread_t threads[nthrds];
  paramsElim *thread_params = (paramsElim *) malloc(nthrds * sizeof(paramsElim));
  uint32 chunkSize;
  uint32 pad;
  uint32 ctr;

  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Naive Gaussian Elimination with pivoting" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();
  for (uint32 i = 0; i < boundary; ++i) {
    A(i,i) %= prime;
#if F4RT_DBG
    std::cout << "!! A(" << i << "," << i << ") " << A(i,i) << std::endl;
    std::cout << "A(" << i << "," << i << ") " << A(i,i) % prime << std::endl;
#endif
    if (A(i,i) == 0) {
      l = i+1;
#if F4RT_DBG
      std::cout << "l1 " << l << std::endl;
#endif
      while (l < m && A(l,i) % prime == 0) {
        l++;
#if F4RT_DBG
        std::cout << "l2 " << l << std::endl;
#endif
      }
      if (l == m) {
        continue;
      } else {
        // swapping
#if F4RT_DBG
        std::cout << "before swapping i = " << i << " with l = " << l << std::endl;
        for (int kk=i; kk < n; ++kk) {
          std::cout << "A(i,"<< kk << ") = " << A(i,kk) << std::endl;
          std::cout << "A(l," << kk << ") = " << A(l,kk) << std::endl;
        }
        std::cout << "n  " << n << std::endl;
        std::cout << "l  " << l << std::endl;
        std::cout << "i  " << i << std::endl;
#endif
        std::vector<mat> tempRow(
            A.entries.begin()+i+(l*n),
            A.entries.begin()+n-1+(l*n)+1);
#if F4RT_DBG
        std::cout << "(n-1)-i " << n-1-i << std::endl;
        std::cout << "sizeof mat " << sizeof(mat) << std::endl;
        std::cout << "sizeof temp " << tempRow.size() << std::endl;
#endif
        for (uint32 it = i; it < n; ++it)
          A.entries[it+l*n] = A.entries[it+i*n];
        for (uint32 it = i; it < n; ++it) {
          A.entries[it+i*n] = tempRow[it-i];
        }
#if F4RT_DBG
        std::cout << "after swapping" << std::endl;
        for (int kk=i; kk < n; ++kk) {
          std::cout << "A(i,"<< kk << ") = " << A(i,kk) << std::endl;
          std::cout << "A(l," << kk << ") = " << A(l,kk) << std::endl;
        }
#endif
        tempRow.clear();
      }
    }
    inv  = negInverseModP(A(i,i), prime);
#if F4RT_DBG
    std::cout << "inv  " << inv << std::endl;
#endif
    chunkSize = (m - i - 1) / nthrds;
    pad       = (m - i - 1) % nthrds;
    ctr = 0;
    for (uint32 l = 0; l < nthrds; ++l) {
      thread_params[l].a      = a_entries;
      thread_params[l].prime  = prime;
      thread_params[l].index  = i;
      if (l < pad) {
        thread_params[l].size   = chunkSize + 1;
        thread_params[l].start  = ctr;
        ctr +=  chunkSize + 1;
      } else {
        thread_params[l].size = chunkSize;
        thread_params[l].start  = ctr;
        ctr +=  chunkSize;
      }
      thread_params[l].n    = n;
      thread_params[l].inv  = inv;
      // real computation
      pthread_create(&threads[l], NULL, elimPTHRD, (void *) &thread_params[l]);
    }

    // join threads back again
    for (uint32 l = 0; l < nthrds; ++l)
      pthread_join(threads[l], NULL);
  }
  free(thread_params);
  //cleanUpModP(A, prime);
  //A.print();
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           pThread 1D" << std::endl;
  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = countGEPFlops(m, n, prime);
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

// cache-oblivious stuff

void elimCoPTHRDBaseModP( mat *M, const uint32 k1, const uint32 i1,
                        const uint32 j1, const uint32 rows, const uint32 cols,
                        uint64 size, uint64 prime, mat *neg_inv_piv, int nthrds) {
  uint64 k;

  for (k = 0; k < size; k++) {
    M[(k1+k)+(k1+k)*cols] %= prime;
    // possibly the negative inverses of the pivots at place (k,k) were already
    // computed in another call. otherwise we need to compute and store it

    if (!neg_inv_piv[k+k1]) {
      if (M[(k1+k)+(k1+k)*cols] != 0) {
        neg_inv_piv[k+k1] = negInverseModP(M[(k1+k)+(k1+k)*cols], prime);
      }
    }
    const mat inv_piv   = neg_inv_piv[k+k1];
    // if the pivots are in the same row part of the matrix as Mmdf then we can
    // always start at the next row (k+1), otherwise we need to start at
    // row 0
    const uint64 istart  = (k1 == i1) ? k+1 : 0;
    //tbb::parallel_for(tbb::blocked_range<uint32>(istart, size, 1),
    //    [&](const tbb::blocked_range<uint32>& r)
    //    {
    for (uint64 i = istart; i < size; i++) {
      //for (uint64 i = r.begin(); i != r.end(); ++i) {
      const mat tmp = (M[k+k1+(i1+i)*cols] * inv_piv) % prime;
      // if the pivots are in the same column part of the matrix as Mmdf then we can
      // always start at the next column (k+1), otherwise we need to start at
      // column 0
      const uint64 jstart  = (k1 == j1) ? k+1 : 0;
      for (uint64 j = jstart; j < size; j++) {
        M[(j1+j)+(i1+i)*cols]  +=  M[(j1+j)+(k1+k)*cols] * tmp;
        M[(j1+j)+(i1+i)*cols]  %=  prime;
      }
    }
    //});
  }
}

void* D1PTHRD(void *p) {
  thrdData *_data       = (thrdData *)p;
  thrdPool *pool        = _data->pool;
  paramsCoElim *params  = _data->params;
  
  if (params->i2 <= params->k1 || params->j2 <= params->k1)
    return 0;

  //printf("[D] Running jobs: %d\n", pool->runningJobs);
  // increase number of active threads
  pthread_mutex_lock(&mutex1);
  pool->runningJobs++;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);
  uint64 size = params->size;

  if (size <= params->blocksize) {
    elimCoPTHRDBaseModP(params->M, params->k1, params->i1, params->j1,
                        params->rows, params->cols, params->size,
                        params->prime, params->neg_inv_piv,
                        params->nthrds);
  } else {
    const uint32 i1 = params->i1;
    const uint32 i2 = params->i2;
    const uint32 j1 = params->j1;
    const uint32 j2 = params->j2;
    const uint32 k1 = params->k1;
    const uint32 k2 = params->k2;

    size = size / 2;

    const uint32 km = (k1+k2) / 2;
    const uint32 im = (i1+i2) / 2;
    const uint32 jm = (j1+j2) / 2;

    pthread_t thread[4];
    paramsCoElim *thread_params = (paramsCoElim *)
                                    malloc(4 * sizeof(paramsCoElim));
    thrdData *data  = (thrdData *) malloc(4 * sizeof(thrdData));
    for (int i = 0; i < 4; ++i) {
      data[i].pool    = pool;
      data[i].params  = &thread_params[i];
    }

    // get not thread-specific parameters -- once for all
    for (int i = 0; i < 4; ++i) {
      thread_params[i].M            = params->M;
      thread_params[i].neg_inv_piv  = params->neg_inv_piv;
      thread_params[i].size         = size;
      thread_params[i].blocksize    = params->blocksize;
      thread_params[i].nthrds       = params->nthrds;
      thread_params[i].prime        = params->prime;
      thread_params[i].rows         = params->rows;
      thread_params[i].cols         = params->cols;
    }

    // X11
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = k1;
    thread_params[0].k2 = km;
    // X12
    thread_params[1].i1 = i1;
    thread_params[1].i2 = im;
    thread_params[1].j1 = jm+1;
    thread_params[1].j2 = j2;
    thread_params[1].k1 = k1;
    thread_params[1].k2 = km;
    // X21
    thread_params[2].i1 = im+1;
    thread_params[2].i2 = i2;
    thread_params[2].j1 = j1;
    thread_params[2].j2 = jm;
    thread_params[2].k1 = k1;
    thread_params[2].k2 = km;
    // X22
    thread_params[3].i1 = im+1;
    thread_params[3].i2 = i2;
    thread_params[3].j1 = jm+1;
    thread_params[3].j2 = j2;
    thread_params[3].k1 = k1;
    thread_params[3].k2 = km;

    // parallel - start
    for (int i = 0; i < 4; ++i)
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 4; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end
    // X11
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = km+1;
    thread_params[0].k2 = k2;
    // X12
    thread_params[1].i1 = i1;
    thread_params[1].i2 = im;
    thread_params[1].j1 = jm+1;
    thread_params[1].j2 = j2;
    thread_params[1].k1 = km+1;
    thread_params[1].k2 = k2;
    // X21
    thread_params[2].i1 = im+1;
    thread_params[2].i2 = i2;
    thread_params[2].j1 = j1;
    thread_params[2].j2 = jm;
    thread_params[2].k1 = km+1;
    thread_params[2].k2 = k2;
    // X22
    thread_params[3].i1 = im+1;
    thread_params[3].i2 = i2;
    thread_params[3].j1 = jm+1;
    thread_params[3].j2 = j2;
    thread_params[3].k1 = km+1;
    thread_params[3].k2 = k2;

    // parallel - start
    for (int i = 0; i < 4; ++i)
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 4; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end
  }
  // decrease number of active threads
  pthread_mutex_lock(&mutex1);
  pool->runningJobs--;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);
  return 0;
}


void* C1PTHRD(void *p) {
  thrdData *_data       = (thrdData *)p;
  thrdPool *pool        = _data->pool;
  paramsCoElim *params  = _data->params;
  if (params->i2 <= params->k1 || params->j2 <= params->k1)
    return 0;

  //printf("[C] Running jobs: %d\n", pool->runningJobs);
  // increase number of active threads
  pthread_mutex_lock(&mutex1);
  pool->runningJobs++;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);
  uint64 size = params->size;

  if (size <= params->blocksize) {
    elimCoPTHRDBaseModP(params->M, params->k1, params->i1, params->j1,
                        params->rows, params->cols, params->size,
                        params->prime, params->neg_inv_piv,
                        params->nthrds);
  } else {
    const uint32 i1 = params->i1;
    const uint32 i2 = params->i2;
    const uint32 j1 = params->j1;
    const uint32 j2 = params->j2;
    const uint32 k1 = params->k1;
    const uint32 k2 = params->k2;

    size = size / 2;

    const uint32 km = (k1+k2) / 2;
    const uint32 im = (i1+i2) / 2;
    const uint32 jm = (j1+j2) / 2;

    pthread_t thread[2];
    paramsCoElim *thread_params = (paramsCoElim *)
                                    malloc(2 * sizeof(paramsCoElim));
    thrdData *data  = (thrdData *) malloc(2 * sizeof(thrdData));
    for (int i = 0; i < 2; ++i) {
      data[i].pool    = pool;
      data[i].params  = &thread_params[i];
    }

    // get not thread-specific parameters -- once for all
    for (int i = 0; i < 2; ++i) {
      thread_params[i].M            = params->M;
      thread_params[i].neg_inv_piv  = params->neg_inv_piv;
      thread_params[i].size         = size;
      thread_params[i].blocksize    = params->blocksize;
      thread_params[i].nthrds       = params->nthrds;
      thread_params[i].prime        = params->prime;
      thread_params[i].rows         = params->rows;
      thread_params[i].cols         = params->cols;
    }

    // X11
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = k1;
    thread_params[0].k2 = km;
    // X21
    thread_params[1].i1 = im+1;
    thread_params[1].i2 = i2;
    thread_params[1].j1 = j1;
    thread_params[1].j2 = jm;
    thread_params[1].k1 = k1;
    thread_params[1].k2 = km;

    // parallel - start
    for (int i = 0; i < 2; ++i)
      pthread_create(&thread[i], NULL, &C1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end

    // X12
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = jm+1;
    thread_params[0].j2 = j2;
    thread_params[0].k1 = k1;
    thread_params[0].k2 = km;
    // X22
    thread_params[1].i1 = im+1;
    thread_params[1].i2 = i2;
    thread_params[1].j1 = jm+1;
    thread_params[1].j2 = j2;
    thread_params[1].k1 = k1;
    thread_params[1].k2 = km;
    // parallel - start
    for (int i = 0; i < 2; ++i)
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end

    // X12
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = jm+1;
    thread_params[0].j2 = j2;
    thread_params[0].k1 = km+1;
    thread_params[0].k2 = k2;
    // X22
    thread_params[1].i1 = im+1;
    thread_params[1].i2 = i2;
    thread_params[1].j1 = jm+1;
    thread_params[1].j2 = j2;
    thread_params[1].k1 = km+1;
    thread_params[1].k2 = k2;
    // parallel - start
    for (int i = 0; i < 2; ++i)
      pthread_create(&thread[i], NULL, &C1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end

    // X11
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = km+1;
    thread_params[0].k2 = k2;
    // X12
    thread_params[1].i1 = im+1;
    thread_params[1].i2 = i2;
    thread_params[1].j1 = j1;
    thread_params[1].j2 = jm;
    thread_params[1].k1 = km+1;
    thread_params[1].k2 = k2;
    // parallel - start
    for (int i = 0; i < 2; ++i)
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end
  }
  // decrease number of active threads
  pthread_mutex_lock(&mutex1);
  pool->runningJobs;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);
  return 0;
}

void* B1PTHRD(void *p) {
  thrdData *_data       = (thrdData *)p;
  thrdPool *pool        = _data->pool;
  paramsCoElim *params  = _data->params;
  if (params->i2 <= params->k1 || params->j2 <= params->k1)
    return 0;

  //printf("[B] Running jobs: %d\n", pool->runningJobs);
  // increase number of active threads
  pthread_mutex_lock(&mutex1);
  pool->runningJobs++;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);
  uint64 size = params->size;

  if (size <= params->blocksize) {
    elimCoPTHRDBaseModP(params->M, params->k1, params->i1, params->j1,
                        params->rows, params->cols, params->size,
                        params->prime, params->neg_inv_piv,
                        params->nthrds);
  } else {
    const uint32 i1 = params->i1;
    const uint32 i2 = params->i2;
    const uint32 j1 = params->j1;
    const uint32 j2 = params->j2;
    const uint32 k1 = params->k1;
    const uint32 k2 = params->k2;

    size = size / 2;

    const uint32 km = (k1+k2) / 2;
    const uint32 im = (i1+i2) / 2;
    const uint32 jm = (j1+j2) / 2;

    pthread_t thread[2];
    paramsCoElim *thread_params = (paramsCoElim *)
                                    malloc(2 * sizeof(paramsCoElim));
    thrdData *data  = (thrdData *) malloc(2 * sizeof(thrdData));
    for (int i = 0; i < 2; ++i) {
      data[i].pool    = pool;
      data[i].params  = &thread_params[i];
    }

    // get not thread-specific parameters -- once for all
    for (int i = 0; i < 2; ++i) {
      thread_params[i].M            = params->M;
      thread_params[i].neg_inv_piv  = params->neg_inv_piv;
      thread_params[i].size         = size;
      thread_params[i].blocksize    = params->blocksize;
      thread_params[i].nthrds       = params->nthrds;
      thread_params[i].prime        = params->prime;
      thread_params[i].rows         = params->rows;
      thread_params[i].cols         = params->cols;
    }

    // X11
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = k1;
    thread_params[0].k2 = km;
    // X12
    thread_params[1].i1 = i1;
    thread_params[1].i2 = im;
    thread_params[1].j1 = jm+1;
    thread_params[1].j2 = j2;
    thread_params[1].k1 = k1;
    thread_params[1].k2 = km;
    // parallel - start
    for (int i = 0; i < 2; ++i)
      pthread_create(&thread[i], NULL, &B1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end

    // X21
    thread_params[0].i1 = im+1;
    thread_params[0].i2 = i2;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = k1;
    thread_params[0].k2 = km;
    // X22
    thread_params[1].i1 = im+1;
    thread_params[1].i2 = i2;
    thread_params[1].j1 = jm+1;
    thread_params[1].j2 = j2;
    thread_params[1].k1 = k1;
    thread_params[1].k2 = km;
    // parallel - start

    for (int i = 0; i < 2; ++i)
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end

    // X21
    thread_params[0].i1 = im+1;
    thread_params[0].i2 = i2;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = km+1;
    thread_params[0].k2 = k2;
    // X22
    thread_params[1].i1 = im+1;
    thread_params[1].i2 = i2;
    thread_params[1].j1 = jm+1;
    thread_params[1].j2 = j2;
    thread_params[1].k1 = km+1;
    thread_params[1].k2 = k2;
    // parallel - start

    for (int i = 0; i < 2; ++i)
      pthread_create(&thread[i], NULL, &B1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end

    // X11
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = km+1;
    thread_params[0].k2 = k2;
    // X12
    thread_params[1].i1 = i1;
    thread_params[1].i2 = im;
    thread_params[1].j1 = jm+1;
    thread_params[1].j2 = j2;
    thread_params[1].k1 = km+1;
    thread_params[1].k2 = k2;
    // parallel - start
    for (int i = 0; i < 2; ++i)
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &data[i]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end
  }
  // decrease number of active threads
  pthread_mutex_lock(&mutex1);
  pool->runningJobs--;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);
  return 0;
}

void* APTHRD(void *p) {
  thrdData *_data       = (thrdData *)p;
  thrdPool *pool        = _data->pool;
  paramsCoElim *params  = _data->params;
  if (params->i2 <= params->k1 || params->j2 <= params->k1)
    return 0;

  //printf("[A] Running jobs: %d\n", pool->runningJobs);
  // increase number of active threads
  pthread_mutex_lock(&mutex1);
  pool->runningJobs++;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);
  uint64 size = params->size;

  if (size <= params->blocksize) {
    elimCoPTHRDBaseModP(params->M, params->k1, params->i1, params->j1,
                        params->rows, params->cols, params->size,
                        params->prime, params->neg_inv_piv,
                        params->nthrds);
  } else {
    const uint32 i1 = params->i1;
    const uint32 i2 = params->i2;
    const uint32 j1 = params->j1;
    const uint32 j2 = params->j2;
    const uint32 k1 = params->k1;
    const uint32 k2 = params->k2;

    size = size / 2;

    const uint32 km = (k1+k2) / 2;
    const uint32 im = (i1+i2) / 2;
    const uint32 jm = (j1+j2) / 2;

    pthread_t thread[2];
    paramsCoElim *thread_params = (paramsCoElim *)
                                    malloc(2 * sizeof(paramsCoElim));
    thrdData *data  = (thrdData *) malloc(2 * sizeof(thrdData));
    for (int i = 0; i < 2; ++i) {
      data[i].pool    = pool;
      data[i].params  = &thread_params[i];
    }

    // get not thread-specific parameters -- once for all
    for (int i = 0; i < 2; ++i) {
      thread_params[i].M            = params->M;
      thread_params[i].neg_inv_piv  = params->neg_inv_piv;
      thread_params[i].size         = size;
      thread_params[i].blocksize    = params->blocksize;
      thread_params[i].nthrds       = params->nthrds;
      thread_params[i].prime        = params->prime;
      thread_params[i].rows         = params->rows;
      thread_params[i].cols         = params->cols;
    }

    // X11
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = j1;
    thread_params[0].j2 = jm;
    thread_params[0].k1 = k1;
    thread_params[0].k2 = km;

    // forward step
    APTHRD((void *) &data[0]);

    // X12
    thread_params[0].i1 = i1;
    thread_params[0].i2 = im;
    thread_params[0].j1 = jm+1;
    thread_params[0].j2 = j2;
    thread_params[0].k1 = k1;
    thread_params[0].k2 = km;
    // X21
    thread_params[1].i1 = im+1;
    thread_params[1].i2 = i2;
    thread_params[1].j1 = j1;
    thread_params[1].j2 = jm;
    thread_params[1].k1 = k1;
    thread_params[1].k2 = km;

    // parallel - start
    if (jobsRunning < maxThreads)
    pthread_create(&thread[0], NULL, &B1PTHRD, (void *) &data[0]);
    pthread_create(&thread[1], NULL, &C1PTHRD, (void *) &data[1]);
    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end

    // get parameters
    thread_params[0].i1 = im+1;
    thread_params[0].i2 = i2;
    thread_params[0].j1 = jm+1;
    thread_params[0].j2 = j2;
    thread_params[0].k1 = k1;
    thread_params[0].k2 = km;

    D1PTHRD((void *) &data[0]);

    // backward step

    thread_params[0].i1 = im+1;
    thread_params[0].i2 = i2;
    thread_params[0].j1 = jm+1;
    thread_params[0].j2 = j2;
    thread_params[0].k1 = km+1;
    thread_params[0].k2 = k2;

    APTHRD((void *) &data[0]);
  }

  // decrease number of active threads
  pthread_mutex_lock(&mutex1);
  pool->runningJobs--;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);
  return 0;
}

void elimCoPTHRDModP(Matrix& M, int nthrds, uint32 blocksize, uint64 prime) {
  if (nthrds <= 0) {
    nthrds  = 1;
  }
  // allocated thread pool
  pthread_t *threads          = (pthread_t *) malloc(nthrds * sizeof(pthread_t));
  thrdPool *pool              = (thrdPool *) malloc(sizeof(thrdPool));
  paramsCoElim *thread_params = (paramsCoElim *)
                                  malloc(sizeof(paramsCoElim));
  thrdData * data             = (thrdData *) malloc(sizeof(thrdData));
  data->pool          = pool;
  data->params        = thread_params;
  pool->maxNumThreads = nthrds;
  pool->runningJobs   = 1;
  printf("[MAIN] Running jobs: %d\n", pool->runningJobs);
  jobsRunning         = 1;
  maxThreads          = nthrds;
  pool->bitmask       = ULONG_MAX;
  // if m > n then only n eliminations are possible
  uint32 m          = M.nRows();
  uint32 n          = M.nCols();
  uint32 boundary   = (m > n) ? n : m;
  mat *a_entries    = M.entries.data();
  mat *neg_inv_piv  = (mat *)calloc(boundary, sizeof(mat));
  a_entries[0]      %=  prime;
  neg_inv_piv[0]    = negInverseModP(a_entries[0], prime);


  printf("Pool: maxNumThreads %d\n", pool->maxNumThreads);
  printf("      runningJobs   %d -- number of bits %d\n", pool->runningJobs, sizeof(pool->runningJobs) * 8);
  printf("      bitmask       %lu -- number of bits %d\n", pool->bitmask, sizeof(pool->bitmask) * 8);
  thread_params->M            = a_entries;
  thread_params->neg_inv_piv  = neg_inv_piv;
  thread_params->blocksize    = blocksize;
  thread_params->nthrds       = nthrds;
  thread_params->rows         = m;
  thread_params->cols         = n;
  thread_params->size         = boundary;
  thread_params->prime        = prime;
  thread_params->i1           = 0;
  thread_params->i2           = m - 1;
  thread_params->j1           = 0;
  thread_params->j2           = n - 1;
  thread_params->k1           = 0;
  thread_params->k2           = boundary - 1;

  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Cache-oblivious Gaussian Elimination without pivoting" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();

  // computation of blocks
  APTHRD((void *) data);
  printf("[MAIN] Running jobs at the end: %d\n", pool->runningJobs);

  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           pThread 1D" << std::endl;
  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = countGEPFlops(m, n, prime);
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
