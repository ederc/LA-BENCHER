#include <matrix.h>

#define __F4RT_DEBUG  0


// multiplies A*B^T and stores it in *this
void multSEQ(Matrix& C, const Matrix& A, const Matrix& B, int blocksize, int impose) {
  // assertion seems strange, but remember that we compute A*B^T
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
  C.resize(l*m);
  std::cout << "Matrix Multiplication" << std::endl;
  timeval start, stop;
  clock_t cStart, cStop;
  float sum = 0;
#if __F4RT_DEBUG
  std::cout << std::endl;
  std::cout << "A => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "B => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "C => " << C.nRows() << "-" << C.nCols() << "-" << C.nEntries() << std::endl;
#endif
  gettimeofday(&start, NULL);
  cStart  = clock();
  if (impose == 1) {
    for (uint32 i = 0; i < l; ++i) {
      for (uint32 j = 0; j < m; ++j) {
        sum = 0;
        for (uint32 k = 0; k < n; k++) {
          // sum += A(i,k) * B(j,k);
          sum += A.entries[k+i*n] * B.entries[k+j*n];
        }
        C.entries[j+i*m]  = (float) (sum);
      }
    }
  } else {
    for (uint32 i = 0; i < l; ++i) {
      for (uint32 j = 0; j < m; ++j) {
        sum = 0;
        for (uint32 k = 0; k < n; k++) {
          // sum += A(i,k) * B(k,j);
          sum += A.entries[k+i*n] * B.entries[j+k*m];
        }
        C.entries[j+i*m]  = (float) (sum);
      }
    }
  }
  gettimeofday(&stop, NULL);
  cStop = clock();
  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  float epsilon = 0.0000000001;
  double realtime = stop.tv_sec - start.tv_sec;
  double cputime  = (cStop - cStart) / CLOCKS_PER_SEC;
  double ratio = cputime/realtime;
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Raw sequential" << std::endl;
  std::cout << "Cache improved:   ";
  if (impose == 1)
    std::cout << "1" << std::endl;
  else
    std::cout << "0" << std::endl;
  std::cout << "# Threads:        1" << std::endl;
  std::cout << "Block size:       " << blocksize << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "Real time:        " << stop.tv_sec - start.tv_sec << " sec" 
    << std::endl;
  std::cout << "CPU time:         " << (cStop - cStart) / CLOCKS_PER_SEC 
    << " sec" << std::    endl;
  if ((int)((cStop - cStart) / CLOCKS_PER_SEC) > epsilon)
    std::cout << "CPU/real time:    " << std::setprecision(4) 
      << ratio << std::endl;
  std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
  std::cout << "GFLOPS/sec:       " << std::setprecision(4) 
    << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std:: endl;
  std::cout << "---------------------------------------------------" << std::endl;
}
