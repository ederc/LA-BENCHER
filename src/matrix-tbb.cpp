#include <matrix.h>

#define BLOCK_SIZE    2
#define __F4RT_DEBUG  0

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
  //
  // do matrix multiplication of submatrices of size in the order of 
  // BLOCK_SIZE
  tbb::parallel_for(tbb::blocked_range<size_t>(0, B.nRows(), BLOCK_SIZE),
      [&](const tbb::blocked_range<size_t>& r)
  {
    for( size_t i=r.begin(); i!=r.end(); ++i )
      for( size_t j=0; j!=B.nRows(); ++j ) {
        uint32 sum = 0;
        for( size_t k=0; k<B.nCols(); ++k )
          sum += A(i,k) * A(j,k);
        C(i,j) = sum;
      }
    
  });

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
