#include <matrix.h>

#define __F4RT_DEBUG  0

// multiplies A*B^T and stores it in *this
void multTBBAuto(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nRows());
  C.resize(A.nRows()*B.nCols());
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
    nthrds  = tbb::task_scheduler_init::default_num_threads();
  }
  tbb::task_scheduler_init init(nthrds);
  // do matrix multiplication of submatrices of size in the order of 
  // blocksize
  gettimeofday(&start, NULL);
  cStart  = clock();
  tbb::parallel_for(tbb::blocked_range<size_t>(0, A.nRows(), blocksize),
      [&](const tbb::blocked_range<size_t>& r)
  {
    for( size_t i=r.begin(); i!=r.end(); ++i )
      for( size_t j=0; j!=B.nCols(); ++j ) {
        float sum = 0;
        for( size_t k=0; k<B.nRows(); ++k )
          sum += A(i,k) * B(k,j);
        C(i,j) = sum;
      }
  });

  gettimeofday(&stop, NULL);
  cStop = clock();

  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "Method:      Intel TBB 1D auto partitioner" << std::endl;
  std::cout << "# Threads:   " << nthrds << std::endl;
  std::cout << "Block size:  " << blocksize << std::endl;
  std::cout << "Real time:   " << stop.tv_sec - start.tv_sec << " sec" << std::endl;
  std::cout << "CPU time:    " << (cStop - cStart) / CLOCKS_PER_SEC << " sec" << std::    endl;
  std::cout << "GFLOPS/sec:  " << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std:: endl;
  std::cout << "----------------------------------------------" << std::endl;
}

// multiplies A*B^T and stores it in *this
void multTBBAffine(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nRows());
  C.resize(A.nRows()*B.nCols());
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
    nthrds  = tbb::task_scheduler_init::default_num_threads();
  }
  tbb::task_scheduler_init init(nthrds);
  tbb::affinity_partitioner ap;
  // do matrix multiplication of submatrices of size in the order of 
  // blocksize
  gettimeofday(&start, NULL);
  cStart  = clock();
  tbb::parallel_for(tbb::blocked_range<size_t>(0, A.nRows(), blocksize),
      [&](const tbb::blocked_range<size_t>& r)
  {
    for( size_t i=r.begin(); i!=r.end(); ++i )
      for( size_t j=0; j!=B.nCols(); ++j ) {
        float sum = 0;
        for( size_t k=0; k<B.nRows(); ++k )
          sum += A(i,k) * B(k,j);
        C(i,j) = sum;
      }
  }, ap);

  gettimeofday(&stop, NULL);
  cStop = clock();

  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "Method:      Intel TBB 1D affinity partitioner" << std::endl;
  std::cout << "# Threads:   " << nthrds << std::endl;
  std::cout << "Block size:  " << blocksize << std::endl;
  std::cout << "Real time:   " << stop.tv_sec - start.tv_sec << " sec" << std::endl;
  std::cout << "CPU time:    " << (cStop - cStart) / CLOCKS_PER_SEC << " sec" << std::    endl;
  std::cout << "GFLOPS/sec:  " << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std:: endl;
  std::cout << "----------------------------------------------" << std::endl;
}

// multiplies A*B^T and stores it in *this
void multTBBSimple(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nRows());
  C.resize(A.nRows()*B.nCols());
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
    nthrds  = tbb::task_scheduler_init::default_num_threads();
  }
  tbb::task_scheduler_init init(nthrds);
  tbb::simple_partitioner sp;
  // do matrix multiplication of submatrices of size in the order of 
  // blocksize
  gettimeofday(&start, NULL);
  cStart  = clock();
  tbb::parallel_for(tbb::blocked_range<size_t>(0, A.nRows(), blocksize),
      [&](const tbb::blocked_range<size_t>& r)
  {
    for( size_t i=r.begin(); i!=r.end(); ++i )
      for( size_t j=0; j!=B.nCols(); ++j ) {
        float sum = 0;
        for( size_t k=0; k<B.nRows(); ++k )
          sum += A(i,k) * B(k,j);
        C(i,j) = sum;
      }
  }, sp);

  gettimeofday(&stop, NULL);
  cStop = clock();

  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "Method:      Intel TBB 1D simple partitioner" << std::endl;
  std::cout << "# Threads:   " << nthrds << std::endl;
  std::cout << "Block size:  " << blocksize << std::endl;
  std::cout << "Real time:   " << stop.tv_sec - start.tv_sec << " sec" << std::endl;
  std::cout << "CPU time:    " << (cStop - cStart) / CLOCKS_PER_SEC << " sec" << std::    endl;
  std::cout << "GFLOPS/sec:  " << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std:: endl;
  std::cout << "----------------------------------------------" << std::endl;
}

// multiplies A*B^T and stores it in *this
void multTBBAuto2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nRows());
  C.resize(A.nRows()*B.nCols());
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
    nthrds  = tbb::task_scheduler_init::default_num_threads();
  }
  tbb::task_scheduler_init init(nthrds);
  // do matrix multiplication of submatrices of size in the order of 
  // blocksize
  gettimeofday(&start, NULL);
  cStart  = clock();
  tbb::parallel_for(tbb::blocked_range2d<size_t>(0, A.nRows(), blocksize, 0, B.nCols(), blocksize),
      [&](const tbb::blocked_range2d<size_t>& r)
  {
    for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i ) {
      for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
        float sum = 0;
        for( size_t k=0; k<B.nRows(); ++k )
          sum += A(i,k) * B(k,j);
        C(i,j) = sum;
      }
    }
  });

  gettimeofday(&stop, NULL);
  cStop = clock();

  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "Method:      Intel TBB 2D auto partitioner" << std::endl;
  std::cout << "# Threads:   " << nthrds << std::endl;
  std::cout << "Block size:  " << blocksize << std::endl;
  std::cout << "Real time:   " << stop.tv_sec - start.tv_sec << " sec" << std::endl;
  std::cout << "CPU time:    " << (cStop - cStart) / CLOCKS_PER_SEC << " sec" << std::    endl;
  std::cout << "GFLOPS/sec:  " << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std:: endl;
  std::cout << "----------------------------------------------" << std::endl;
}

// multiplies A*B^T and stores it in *this
void multTBBAffine2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nRows());
  C.resize(A.nRows()*B.nCols());
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
    nthrds  = tbb::task_scheduler_init::default_num_threads();
  }
  tbb::task_scheduler_init init(nthrds);
  tbb::affinity_partitioner ap;
  // do matrix multiplication of submatrices of size in the order of 
  // blocksize
  gettimeofday(&start, NULL);
  cStart  = clock();
  tbb::parallel_for(tbb::blocked_range2d<size_t>(0, A.nRows(), blocksize, 0, B.nCols(), blocksize),
      [&](const tbb::blocked_range2d<size_t>& r)
  {
    for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i ) {
      for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
        float sum = 0;
        for( size_t k=0; k<B.nRows(); ++k )
          sum += A(i,k) * B(k,j);
        C(i,j) = sum;
      }
    }
  }, ap);

  gettimeofday(&stop, NULL);
  cStop = clock();

  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "Method:      Intel TBB 2D affinity partitioner" << std::endl;
  std::cout << "# Threads:   " << nthrds << std::endl;
  std::cout << "Block size:  " << blocksize << std::endl;
  std::cout << "Real time:   " << stop.tv_sec - start.tv_sec << " sec" << std::endl;
  std::cout << "CPU time:    " << (cStop - cStart) / CLOCKS_PER_SEC << " sec" << std::    endl;
  std::cout << "GFLOPS/sec:  " << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std:: endl;
  std::cout << "----------------------------------------------" << std::endl;
}

// multiplies A*B^T and stores it in *this
void multTBBSimple2d(Matrix& C, const Matrix& A, const Matrix& B, int nthrds, int blocksize) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nRows());
  C.resize(A.nRows()*B.nCols());
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
    nthrds  = tbb::task_scheduler_init::default_num_threads();
  }
  tbb::task_scheduler_init init(nthrds);
  tbb::simple_partitioner sp;
  // do matrix multiplication of submatrices of size in the order of 
  // blocksize
  gettimeofday(&start, NULL);
  cStart  = clock();
  tbb::parallel_for(tbb::blocked_range2d<size_t>(0, A.nRows(), blocksize, 0, B.nCols(), blocksize),
      [&](const tbb::blocked_range2d<size_t>& r)
  {
    for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i ) {
      for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
        float sum = 0;
        for( size_t k=0; k<B.nRows(); ++k )
          sum += A(i,k) * B(k,j);
        C(i,j) = sum;
      }
    }
  }, sp);

  gettimeofday(&stop, NULL);
  cStop = clock();

  // compute FLOPS:
  // assume addition and multiplication in the mult kernel are 2 operations
  // done A.nRows() * B.nRows() * B.nCols()
  double flops = 2 * A.nRows() * B.nRows() * B.nCols();
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "Method:      Intel TBB 2D simple partitioner" << std::endl;
  std::cout << "# Threads:   " << nthrds << std::endl;
  std::cout << "Block size:  " << blocksize << std::endl;
  std::cout << "Real time:   " << stop.tv_sec - start.tv_sec << " sec" << std::endl;
  std::cout << "CPU time:    " << (cStop - cStart) / CLOCKS_PER_SEC << " sec" << std::    endl;
  std::cout << "GFLOPS/sec:  " << flops / (1000000000 * (stop.tv_sec - start.tv_sec)) << std:: endl;
  std::cout << "----------------------------------------------" << std::endl;
}
