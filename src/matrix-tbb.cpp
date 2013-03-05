#include <matrix-tbb.h>

// multiplies A*B^T and stores it in *this
void multTBB(Matrix& C, const Matrix& A, const Matrix& B) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nCols());
  C.resize(A.nRows()*B.nRows());
  std::cout << " --- Start Matrix Multiplication --- " << std::endl;
  timeval start, stop;
  clock_t cStart, cStop;
  gettimeofday(&start, NULL);
  cStart  = clock();
#if __F4RT_DEBUG
  std::cout << "A => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "B => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "C => " << C.nRows() << "-" << C.nCols() << "-" << C.nEntries() << std::endl;
#endif
  std::vector<uint32> sum;
  const int padding = __F4RT_CPU_CACHE_LINE / sizeof(uint32);
  std::cout << "padding " << padding << std::endl;
  sum.resize(padding * A.nRows() * B.nRows());
  
  // do matrix multiplication of submatrices of size in the order of 32x32
  parallel_for  ( blocked_range2d<size_t>(0, B.nRows(), 32, 0 , B.nRows(), 32), 
                  MatrixMultiply2D(C,A,B)
                );
  gettimeofday(&stop, NULL);
  cStop = clock();

  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  std::cout << " --- End Matrix Multiplication --- " << std::endl;
  std::cout << "REAL TIME for Matrix Multiplication: " << stop.tv_sec - start.tv_sec << " seconds." << std::endl;
  std::cout << "CPU  TIME for Matrix Multiplication: " << (cStop - cStart) / CLOCKS_PER_SEC << " seconds." << std::endl;
  std::cout << "GFLOPS:                              " << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std::endl;
}
