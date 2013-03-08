#include <matrix.h>

#define __F4RT_DEBUG  0


// multiplies A*B^T and stores it in *this
void multSEQ(Matrix& C, const Matrix& A, const Matrix& B, int blocksize) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nRows());
  C.resize(A.nRows()*B.nCols());
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
  for (uint32 i = 0; i < A.nRows(); ++i) {
    for (uint32 j = 0; j < B.nCols(); ++j) {
      float sum = 0;
      for (uint32 k = 0; k < B.nRows(); k++) {
        sum += A(i,k) * B(k,j);
      }
      C(i,j)  = sum;
    }
  }
  gettimeofday(&stop, NULL);
  cStop = clock();
  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "Method:      Sequential" << std::endl;
  std::cout << "# Threads:   1" << std::endl;
  std::cout << "Block size:  " << blocksize << std::endl;
  std::cout << "Real time:   " << stop.tv_sec - start.tv_sec << " sec" << std::endl;
  std::cout << "CPU time:    " << (cStop - cStart) / CLOCKS_PER_SEC << " sec" << std::    endl;
  std::cout << "GFLOPS/sec:  " << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std:: endl;
  std::cout << "----------------------------------------------" << std::endl;
}
