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

void elimNaivePTHRDModP1d(Matrix& A, int nthrds, uint32 blocksize, uint64 prime) {
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
    for (int l = 0; l < nthrds; ++l) {
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
    for (int l = 0; l < nthrds; ++l)
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

void elimNaivePTHRDModP1dPivot(Matrix& A, int nthrds, uint32 blocksize, uint64 prime) {
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
    for (int l = 0; l < nthrds; ++l) {
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
    for (int l = 0; l < nthrds; ++l)
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
  paramsCoElim *_p  = (paramsCoElim *)p;
  // increase number of active threads
  pthread_mutex_lock(&mutex1);
  counter++;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);

  if (_p->i2 <= _p->k1 || _p->j2 <= _p->k1)
    return 0;

  uint64 size = _p->size;

  if (size <= _p->blocksize) {
    elimCoPTHRDBaseModP(_p->M, _p->k1, _p->i1, _p->j1, _p->rows, _p->cols, 
                        _p->size, _p->prime, _p->neg_inv_piv, _p->nthrds);
  } else {
    const uint32 i1 = _p->i1; 
    const uint32 i2 = _p->i2; 
    const uint32 j1 = _p->j1; 
    const uint32 j2 = _p->j2; 
    const uint32 k1 = _p->k1; 
    const uint32 k2 = _p->k2;

    size = size / 2;

    const uint32 km = (k1+k2) / 2;
    const uint32 im = (i1+i2) / 2;
    const uint32 jm = (j1+j2) / 2;

    pthread_t thread[4];
    paramsCoElim *thread_params = (paramsCoElim *) 
                                    malloc(4 * sizeof(paramsCoElim));


    // get parameters
    for (int i = 0; i < 4; ++i) {
      thread_params[i].M            = _p->M;
      thread_params[i].neg_inv_piv  = _p->neg_inv_piv;
      thread_params[i].blocksize    = _p->blocksize;
      thread_params[i].size         = size;
      thread_params[i].nthrds       = _p->nthrds;
      thread_params[i].prime        = _p->prime;
      thread_params[i].rows         = _p->rows;
      thread_params[i].cols         = _p->cols;
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
    for (int i = 0; i < 4; ++i) {
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &thread_params[i]);
    }
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
    for (int i = 0; i < 4; ++i) {
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &thread_params[i]);
    }
    for (int i = 0; i < 4; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end
  }

  // decrease number of active threads
  pthread_mutex_lock(&mutex1);
  counter--;
  //printf("Counter value: %d\n",counter);
  pthread_mutex_unlock(&mutex1);

  return 0;
}

void C1PTHRD(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize) {
  if (i2 <= k1 || j2 <= k1)
    return;

  if (size <= blocksize) {
    elimCoPTHRDBaseModP(M, k1, i1, j1, rows, cols, size, prime,
                      neg_inv_piv, nthrds);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // parallel - start
    // X11
    C1PTHRD( M, k1, km, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
    // X21
    C1PTHRD( M, k1, km, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
    // parallel - end

    pthread_t thread[2];
    paramsCoElim *thread_params = (paramsCoElim *) 
                                    malloc(2 * sizeof(paramsCoElim));


    // get parameters
    for (int i = 0; i < 2; ++i) {
      thread_params[i].M            = M;
      thread_params[i].neg_inv_piv  = neg_inv_piv;
      thread_params[i].size         = size;
      thread_params[i].blocksize    = blocksize;
      thread_params[i].nthrds       = nthrds;
      thread_params[i].prime        = prime;
      thread_params[i].rows         = rows;
      thread_params[i].cols         = cols;
    }

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
    for (int i = 0; i < 2; ++i) {
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &thread_params[i]);
    }

    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end

    // parallel - start
    // X12
    C1PTHRD( M, km+1, k2, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
    // X22
    C1PTHRD( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
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
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &thread_params[i]);

    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end
  }
}

void B1PTHRD(mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize) {
  if (i2 <= k1 || j2 <= k1)
    return;

  if (size <= blocksize) {
    elimCoPTHRDBaseModP(M, k1, i1, j1, rows, cols, size, prime,
                      neg_inv_piv, nthrds);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // parallel - start
    // X11
    B1PTHRD( M, k1, km, i1, im, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
    // X11
    B1PTHRD( M, k1, km, i1, im, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
    // parallel - end

    pthread_t thread[2];
    paramsCoElim *thread_params = (paramsCoElim *) 
                                    malloc(2 * sizeof(paramsCoElim));


    // get parameters
    for (int i = 0; i < 2; ++i) {
      thread_params[i].M            = M;
      thread_params[i].neg_inv_piv  = neg_inv_piv;
      thread_params[i].size         = size;
      thread_params[i].blocksize    = blocksize;
      thread_params[i].nthrds       = nthrds;
      thread_params[i].prime        = prime;
      thread_params[i].rows         = rows;
      thread_params[i].cols         = cols;
    }

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
    
    for (int i = 0; i < 2; ++i) {
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &thread_params[i]);
    }

    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end
    // parallel - start
    // X21
    B1PTHRD( M, km+1, k2, im+1, i2, j1, jm, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
    // X22
    B1PTHRD( M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
        prime, neg_inv_piv, nthrds, blocksize);
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
      pthread_create(&thread[i], NULL, &D1PTHRD, (void *) &thread_params[i]);

    for (int i = 0; i < 2; ++i)
      pthread_join(thread[i], NULL);
    // parallel - end
  }
}

void APTHRD( mat *M, const uint32 k1, const uint32 k2,
        const uint32 i1, const uint32 i2,
		    const uint32 j1, const uint32 j2,
		    const uint32 rows, const uint32 cols,
        uint64 size, uint64 prime, mat *neg_inv_piv,
        int nthrds, uint32 blocksize) {
  if (i2 <= k1 || j2 <= k1)
    return;

  if (size <= blocksize) {
    elimCoPTHRDBaseModP(M, k1, i1, j1, rows, cols, size, prime,
                      neg_inv_piv, nthrds);
  } else {
    size = size / 2;

    uint32 km = (k1+k2) / 2 ;
    uint32 im = (i1+i2) / 2;
    uint32 jm = (j1+j2) / 2;

    // forward step
    APTHRD(M, k1, km, i1, im, j1, jm, rows, cols, size,
      prime, neg_inv_piv, nthrds, blocksize);
/*
    PTHRD_spawn(0, 14, APTHRD,
        PTHRD_MODE_RW, PTHRD_TYPE_INT, rows*cols, M,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, k1,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, km,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, i1,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, im,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, j1,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, jm,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, rows,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, cols,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, size,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, prime,
        //PTHRD_MODE_RW, PTHRD_TYPE_INT, sizeof(neg_inv_piv)/sizeof(neg_inv_piv[0]), neg_inv_piv,
        PTHRD_MODE_RW, PTHRD_TYPE_INT, 4, neg_inv_piv,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, nthrds,
        PTHRD_MODE_V, PTHRD_TYPE_INT, 1, blocksize);
*/
    // parallel - start
    B1PTHRD(M, k1, km, i1, im, jm+1, j2, rows, cols, size,
          prime, neg_inv_piv, nthrds, blocksize);
    C1PTHRD(M, k1, km, im+1, i2, j1, jm, rows, cols, size,
          prime, neg_inv_piv, nthrds, blocksize);
    // parallel - end

    pthread_t thread;
    paramsCoElim *thread_params = (paramsCoElim *) 
                                    malloc(sizeof(paramsCoElim));


    // get parameters
    thread_params->M            = M;
    thread_params->neg_inv_piv  = neg_inv_piv;
    thread_params->size         = size;
    thread_params->blocksize    = blocksize;
    thread_params->nthrds       = nthrds;
    thread_params->prime        = prime;
    thread_params->rows         = rows;
    thread_params->cols         = cols;
    thread_params->i1           = im+1;
    thread_params->i2           = i2;
    thread_params->j1           = jm+1;
    thread_params->j2           = j2;
    thread_params->k1           = k1;
    thread_params->k2           = km;

    D1PTHRD((void *) thread_params);
    //pthread_create(&thread, NULL, &D1PTHRD, (void *) &thread_params);
    //pthread_join(thread, NULL);

    //D1PTHRD( M, k1, km, im+1, i2, jm+1, j2, rows, cols, size,
    //    prime, neg_inv_piv, nthrds, blocksize);

    // backward step
    
    APTHRD(M, km+1, k2, im+1, i2, jm+1, j2, rows, cols, size,
      prime, neg_inv_piv, nthrds, blocksize);
  }
}

void elimCoPTHRDModP(Matrix& M, int nthrds, uint32 blocksize, uint64 prime) {
  uint32 m          = M.nRows();
  uint32 n          = M.nCols();
  // if m > n then only n eliminations are possible
  uint32 boundary   = (m > n) ? n : m;
  mat *a_entries    = M.entries.data();
  mat *neg_inv_piv  =   (mat *)calloc(boundary, sizeof(mat));
  a_entries[0]      %=  prime;
  neg_inv_piv[0]    =   negInverseModP(a_entries[0], prime);
  timeval start, stop;
  clock_t cStart, cStop;
  std::cout << "Cache-oblivious Gaussian Elimination without pivoting" << std::endl;
  gettimeofday(&start, NULL);
  cStart  = clock();

  // computation of blocks
  APTHRD( a_entries, 0, boundary-1, 0, m-1, 0, n-1, m, n,
        boundary, prime, neg_inv_piv, nthrds, blocksize);

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
