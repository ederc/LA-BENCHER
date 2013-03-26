#include <f4rt-omp.h>

#define __F4RT_DEBUG  0

#ifdef __F4RT_HAVE_OPENMP
// multiplies A*B^T and stores it in *this
void multOMP1dInner(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose) {
  uint32 l, m, n;
  int thrdCounter = nthrds;
  if (impose == 1) {
    l = A.nRows();
    m = B.nRows();
    n = B.nCols();
  } else {
    l = A.nRows();
    m = B.nCols();
    n = B.nRows();
  }
  C.resize(l*m);
  std::cout << "Matrix Multiplication" << std::endl;
  timeval start, stop;
  clock_t cStart, cStop;
#if __F4RT_DEBUG
  std::cout << std::endl;
  std::cout << "A => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "B => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "C => " << C.nRows() << "-" << C.nCols() << "-" << C.nEntries() << std::endl;
#endif
  gettimeofday(&start, NULL);
  cStart  = clock();
  if (nthrds > 0)
    omp_set_num_threads(nthrds);
if (impose == 1) {
#pragma omp parallel shared(A,B,C)
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
  for (uint32 i = 0; i < l; ++i) {
#pragma omp for schedule(guided,blocksize) collapse(1)
    for (uint32 j = 0; j < m; ++j) {
      float sum = 0;
      for (uint32 k = 0; k < n; k++) {
        // sum += A(i,k) * B(j,k);
        sum += A.entries[k+i*n] * B.entries[k+j*n];
      }
      C.entries[j+i*m]  = (float) (sum);
    }
  }
}
} else {
#pragma omp parallel shared(A,B,C)
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
#pragma omp for schedule(guided,blocksize) collapse(2)
  for (uint32 i = 0; i < l; ++i) {
    for (uint32 j = 0; j < m; ++j) {
      float sum = 0;
      for (uint32 k = 0; k < n; k++) {
        // sum += A(i,k) * B(k,j);
        sum += A.entries[k+i*n] * B.entries[j+k*m];
      }
      C.entries[j+i*m]  = (float) (sum);
    }
  }
}
}
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Open MP collapse(1) inner loop" << std::endl;
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
  std::cout << "# Threads:        " << thrdCounter << std::endl;
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

// multiplies A*B^T and stores it in *this
void multOMP1dOuter(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose) {
  uint32 l, m, n;
  int thrdCounter = nthrds;
  if (impose == 1) {
    l = A.nRows();
    m = B.nRows();
    n = B.nCols();
  } else {
    l = A.nRows();
    m = B.nCols();
    n = B.nRows();
  }
  C.resize(l*m);
  std::cout << "Matrix Multiplication" << std::endl;
  timeval start, stop;
  clock_t cStart, cStop;
#if __F4RT_DEBUG
  std::cout << std::endl;
  std::cout << "A => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "B => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "C => " << C.nRows() << "-" << C.nCols() << "-" << C.nEntries() << std::endl;
#endif
  gettimeofday(&start, NULL);
  cStart  = clock();
  if (nthrds > 0)
    omp_set_num_threads(nthrds);
if (impose == 1) {
#pragma omp parallel shared(A,B,C)
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
#pragma omp for schedule(guided,blocksize) collapse(1)
  for (uint32 i = 0; i < l; ++i) {
    for (uint32 j = 0; j < m; ++j) {
      float sum = 0;
      for (uint32 k = 0; k < n; k++) {
        // sum += A(i,k) * B(j,k);
        sum += A.entries[k+i*n] * B.entries[k+j*n];
      }
      C.entries[j+i*m]  = (float) (sum);
    }
  }
}
} else {
#pragma omp parallel shared(A,B,C)
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
#pragma omp for schedule(guided,blocksize) collapse(2)
  for (uint32 i = 0; i < l; ++i) {
    for (uint32 j = 0; j < m; ++j) {
      float sum = 0;
      for (uint32 k = 0; k < n; k++) {
        // sum += A(i,k) * B(k,j);
        sum += A.entries[k+i*n] * B.entries[j+k*m];
      }
      C.entries[j+i*m]  = (float) (sum);
    }
  }
}
}
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Open MP collapse(1) outer loop" << std::endl;
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
  std::cout << "# Threads:        " << thrdCounter << std::endl;
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

// multiplies A*B^T and stores it in *this
void multOMP2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize, int impose) {
  uint32 l, m, n;
  int thrdCounter = nthrds;
  if (impose == 1) {
    l = A.nRows();
    m = B.nRows();
    n = B.nCols();
  } else {
    l = A.nRows();
    m = B.nCols();
    n = B.nRows();
  }
  C.resize(l*m);
  std::cout << "Matrix Multiplication" << std::endl;
  timeval start, stop;
  clock_t cStart, cStop;
#if __F4RT_DEBUG
  std::cout << std::endl;
  std::cout << "A => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "B => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "C => " << C.nRows() << "-" << C.nCols() << "-" << C.nEntries() << std::endl;
#endif
  gettimeofday(&start, NULL);
  cStart  = clock();
  if (nthrds > 0)
    omp_set_num_threads(nthrds);
if (impose == 1) {
#pragma omp parallel shared(A,B,C)
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
#pragma omp for schedule(guided,blocksize) collapse(2)
  for (uint32 i = 0; i < l; ++i) {
    for (uint32 j = 0; j < m; ++j) {
      float sum = 0;
      for (uint32 k = 0; k < n; k++) {
        // sum += A(i,k) * B(j,k);
        sum += A.entries[k+i*n] * B.entries[k+j*n];
      }
      C.entries[j+i*m]  = (float) (sum);
    }
  }
}
} else {
#pragma omp parallel shared(A,B,C)
{
  #pragma omp master
  {
    thrdCounter = omp_get_num_threads();
  }
#pragma omp for schedule(guided,blocksize) collapse(2)
  for (uint32 i = 0; i < l; ++i) {
    for (uint32 j = 0; j < m; ++j) {
      float sum = 0;
      for (uint32 k = 0; k < n; k++) {
        // sum += A(i,k) * B(k,j);
        sum += A.entries[k+i*n] * B.entries[j+k*m];
      }
      C.entries[j+i*m]  = (float) (sum);
    }
  }
}
}
  gettimeofday(&stop, NULL);
  cStop = clock();
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Open MP collapse(2)" << std::endl;
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
  std::cout << "# Threads:        " << thrdCounter << std::endl;
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
