#include <matrix.h>

#define __F4RT_DEBUG  0

#ifdef __F4RT_HAVE_INTEL_TBB
// multiplies A*B^T and stores it in *this
void multTBBAuto( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                  int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, l, blocksize),
        [&](const tbb::blocked_range<size_t>& r)
        {
          for( size_t i=r.begin(); i!=r.end(); ++i )
            for( size_t j=0; j!=m; ++j ) {
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              C.entries[j+i*m]  = (float) (sum);
            }
        });
  } else {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, l, blocksize),
        [&](const tbb::blocked_range<size_t>& r)
        {
          for( size_t i=r.begin(); i!=r.end(); ++i )
            for( size_t j=0; j!=m; ++j ) {
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              C.entries[j+i*m]  = (float) (sum);
            }
        });
  }

  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 1D auto partitioner" << std::endl;
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

// multiplies A*B^T and stores it in *this
void multTBBAffine( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                    int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, l, blocksize),
        [&](const tbb::blocked_range<size_t>& r)
        {
          for( size_t i=r.begin(); i!=r.end(); ++i )
            for( size_t j=0; j!=m; ++j ) {
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, ap);
  } else {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, l, blocksize),
        [&](const tbb::blocked_range<size_t>& r)
        {
          for( size_t i=r.begin(); i!=r.end(); ++i )
            for( size_t j=0; j!=m; ++j ) {
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, ap);
  }
  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 1D affinity partitioner" << std::endl;
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

// multiplies A*B^T and stores it in *this
void multTBBSimple( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                    int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, l, blocksize),
        [&](const tbb::blocked_range<size_t>& r)
        {
          for( size_t i=r.begin(); i!=r.end(); ++i )
            for( size_t j=0; j!=m; ++j ) {
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, sp);
  } else {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, l, blocksize),
        [&](const tbb::blocked_range<size_t>& r)
        {
          for( size_t i=r.begin(); i!=r.end(); ++i )
            for( size_t j=0; j!=m; ++j ) {
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, sp);
  }

  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 1D simple partitioner" << std::endl;
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

// multiplies A*B^T and stores it in *this
void multTBBAuto2d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                    int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range2d<size_t>(0, l, blocksize, 0,
                      m, blocksize),
        [&](const tbb::blocked_range2d<size_t>& r)
        {
          for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i )
            for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              C.entries[j+i*m]  = (float) (sum);
            }
        });
  } else {
    tbb::parallel_for(tbb::blocked_range2d<size_t>(0, l, blocksize, 0,
                      m, blocksize),
        [&](const tbb::blocked_range2d<size_t>& r)
        {
          for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i )
            for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              C.entries[j+i*m]  = (float) (sum);
            }
        });
  }

  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 2D auto partitioner" << std::endl;
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

// multiplies A*B^T and stores it in *this
void multTBBAffine2d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                      int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range2d<size_t>(0, l, blocksize, 0,
                      m, blocksize),
        [&](const tbb::blocked_range2d<size_t>& r)
        {
          for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i )
            for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, ap);
  } else {
    tbb::parallel_for(tbb::blocked_range2d<size_t>(0, l, blocksize, 0,
                      m, blocksize),
        [&](const tbb::blocked_range2d<size_t>& r)
        {
          for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i )
            for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
              float sum = 0;
              for( size_t k=0; k<n; ++k )
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, ap);
  }
  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 2D affinity partitioner" << std::endl;
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

// multiplies A*B^T and stores it in *this
void multTBBSimple2d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds,
                      int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range2d<size_t>(0, l, blocksize, 0,
                      m, blocksize),
        [&](const tbb::blocked_range2d<size_t>& r)
        {
          for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i )
            for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
              float sum = 0;
              for( size_t k=0; k<n; ++k ) {
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              }
              //std::cout << j+i*m << sum << std::endl;
              C.entries[j+i*m]  = (float) (sum);
            }
        }, sp);
  } else {
    tbb::parallel_for(tbb::blocked_range2d<size_t>(0, l, blocksize, 0,
                      m, blocksize),
        [&](const tbb::blocked_range2d<size_t>& r)
        {
          for( size_t i=r.rows().begin(); i!=r.rows().end(); ++i )
            for( size_t j=r.cols().begin(); j!=r.cols().end(); ++j ){
              float sum = 0;
              for( size_t k=0; k<n; ++k ) {
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              }
              C.entries[j+i*m]  = (float) (sum);
            }
        }, sp);
  }
  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 2D simple partitioner" << std::endl;
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

// multiplies A*B^T and stores it in *this
void multTBBAuto3d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                    int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range3d<size_t>(0, l, blocksize, 0,
                      m, blocksize, 0, n, blocksize),
        [&](const tbb::blocked_range3d<size_t>& r)
        {
          for( size_t i=r.pages().begin(); i!=r.pages().end(); ++i )
            for( size_t j=r.rows().begin(); j!=r.rows().end(); ++j ){
              float sum = 0;
              for( size_t k=r.cols().begin(); k!=r.cols().end(); ++k )
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              C.entries[j+i*m]  = (float) (sum);
            }
        });
  } else {
    tbb::parallel_for(tbb::blocked_range3d<size_t>(0, l, blocksize, 0,
                      m, blocksize, 0, n, blocksize),
        [&](const tbb::blocked_range3d<size_t>& r)
        {
          for( size_t i=r.pages().begin(); i!=r.pages().end(); ++i )
            for( size_t j=r.rows().begin(); j!=r.rows().end(); ++j ){
              float sum = 0;
              for( size_t k=r.cols().begin(); k!=r.cols().end(); ++k )
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              C.entries[j+i*m]  = (float) (sum);
            }
        });
  }

  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 3D auto partitioner" << std::endl;
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

// multiplies A*B^T and stores it in *this
void multTBBAffine3d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds, 
                      int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range3d<size_t>(0, l, blocksize, 0,
                      m, blocksize, 0, n, blocksize),
        [&](const tbb::blocked_range3d<size_t>& r)
        {
          for( size_t i=r.pages().begin(); i!=r.pages().end(); ++i )
            for( size_t j=r.rows().begin(); j!=r.rows().end(); ++j ){
              float sum = 0;
              for( size_t k=r.cols().begin(); k!=r.cols().end(); ++k )
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, ap);
  } else {
    tbb::parallel_for(tbb::blocked_range3d<size_t>(0, l, blocksize, 0,
                      m, blocksize, 0, n, blocksize),
        [&](const tbb::blocked_range3d<size_t>& r)
        {
          for( size_t i=r.pages().begin(); i!=r.pages().end(); ++i )
            for( size_t j=r.rows().begin(); j!=r.rows().end(); ++j ){
              float sum = 0;
              for( size_t k=r.cols().begin(); k!=r.cols().end(); ++k )
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, ap);
  }
  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 3D affinity partitioner" << std::endl;
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

// multiplies A*B^T and stores it in *this
void multTBBSimple3d( Matrix& C, const Matrix& A, const Matrix& B, int nthrds,
                      int blocksize, int impose) {
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
  if (impose == 1) {
    tbb::parallel_for(tbb::blocked_range3d<size_t>(0, l, blocksize, 0,
                      m, blocksize, 0, n, blocksize),
        [&](const tbb::blocked_range3d<size_t>& r)
        {
          for( size_t i=r.pages().begin(); i!=r.pages().end(); ++i )
            for( size_t j=r.rows().begin(); j!=r.rows().end(); ++j ){
              float sum = 0;
              for( size_t k=r.cols().begin(); k!=r.cols().end(); ++k )
                sum += A.entries[k+i*n] * B.entries[k+j*n];
              //std::cout << j+i*m << sum << std::endl;
              C.entries[j+i*m]  = (float) (sum);
            }
        }, sp);
  } else {
    tbb::parallel_for(tbb::blocked_range3d<size_t>(0, l, blocksize, 0,
                      m, blocksize, 0, n, blocksize),
        [&](const tbb::blocked_range3d<size_t>& r)
        {
          for( size_t i=r.pages().begin(); i!=r.pages().end(); ++i )
            for( size_t j=r.rows().begin(); j!=r.rows().end(); ++j ){
              float sum = 0;
              for( size_t k=r.cols().begin(); k!=r.cols().end(); ++k )
                sum += A.entries[k+i*n] * B.entries[j+k*m];
              C.entries[j+i*m]  = (float) (sum);
            }
        }, sp);
  }
  gettimeofday(&stop, NULL);
  cStop = clock();

  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "Method:           Intel TBB 3D simple partitioner" << std::endl;
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
