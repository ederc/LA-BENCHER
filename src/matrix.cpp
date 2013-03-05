#include "matrix.h"

#define __F4RT_DEBUG  0
namespace {
  template<class T>
  T readOne(FILE* file) {
    T t;
    if (fread(&t, sizeof(T), 1, file) != 1)
      std::cerr << "Error while reading file." << std::endl;
    return t;
  }

  template<class T>
  void readMany(FILE* file, size_t count, std::vector<T>& v) {
    //size_t const origSize = v.size();
    //v.resize(origSize+count);
    if (fread(v.data(), sizeof(T), v.size(), file) != v.size())
      std::cerr << "Error while reading file." << std::endl;
  }

  template<class T>
  void writeOne(const T& t, FILE* file) {
    if (fwrite(&t, sizeof(T), 1, file) != 1)
      std::cerr << "Error while writing to file." << std::endl;
  }

  template<class T>
  void writeMany(const std::vector<T>& v, FILE* file) {
    if (v.empty())
      return;
    if (fwrite(v.data(), sizeof(T), v.size(), file) != v.size())
      std::cerr << "Error while writing to file." << std::endl;
  }

  uint32 getRandom() {
    return std::rand()%65536;
  }

  void check(const Matrix& A, const Matrix& B) {
    std::cout << "Checking generated matrix with data written to file." << std::endl;
    for (size_t i = 0; i < A.nRows(); ++i)
      for (size_t j = 0; j < A.nCols(); ++j) 
        if (A(i,j) != B(i,j))
        std::cerr << "Error: Matrices do not coincide at position (" << i << "," << j << ")" << std::endl;
  }
}

void Matrix::read(FILE* file) {
  const auto rowCount   = readOne<uint32>(file);
  const auto colCount   = readOne<uint32>(file);
  const auto entryCount = rowCount * colCount;
  m = rowCount;
  n = colCount;
  l = entryCount;
  entries.resize(entryCount);
  readMany(file, entryCount, entries);
}

void Matrix::write(FILE* file) {
  writeOne(static_cast<uint32>(nRows()), file);
  writeOne(static_cast<uint32>(nCols()), file);
  writeMany(entries, file);
}

void Matrix::print() {
  std::cout << "#rows" << std::setw(10) << m << std::endl;
  std::cout << "#cols" << std::setw(10) << n << std::endl;
  if (m*n > 1000) {
    int pr;
    std::cout << "NOTE: The matrix consists of " << m*n << " entries." << std::endl;
    std::cout << "      Do you really want to print it? (1=yes, else=no) ";
    std::cin >> pr;
    if (pr != 1)
      return;
  }
  uint32 i,j;
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n-1; ++j) {
      std::cout  << std::setw(5) << entries[i+m*j] << " | ";
    }
    std::cout  << std::setw(5) << entries[i+m*j];
    std::cout << std::endl;
  }
}

// multiplies A*B^T and stores it in *this
void Matrix::mult(const Matrix& A, const Matrix& B) {
  // assertion seems strange, but remember that we compute A*B^T
  assert (A.nCols() == B.nCols());
  m = A.nRows();
  n = B.nRows();
  l = nRows()*nCols();
  entries.resize(l);
  std::cout << " --- Start Matrix Multiplication --- " << std::endl;
  timeval start, stop;
  clock_t cStart, cStop;
  gettimeofday(&start, NULL);
  cStart  = clock();
#if __F4RT_DEBUG
  std::cout << "A => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "B => " << A.nRows() << "-" << A.nCols() << "-" << A.nEntries() << std::endl;
  std::cout << "C => " << nRows() << "-" << nCols() << "-" << nEntries() << std::endl;
#endif
  std::vector<uint32> sum;
  const int padding = __F4RT_CPU_CACHE_LINE / sizeof(uint32);
  std::cout << "padding " << padding << std::endl;
  sum.resize(padding*m*n);
#pragma omp parallel shared(sum)
{
#pragma omp for
  for (uint32 i = 0; i < A.nRows(); ++i) {
    for (uint32 j = 0; j < B.nRows(); ++j) {
#if __F4RT_DEBUG
      std::cout << i << "--" << j << std::endl;
#endif
      for (uint32 k = 0; k < B.nCols(); k++) {
#if __F4RT_DEBUG
        std::cout << "A: " << A(i,k) << " B: " << B(j,k) << std::endl;
#endif
        sum[padding*(i + (j*m))]  += A(i,k) * B(j,k);
      }
      (*this)(i,j)  = sum[padding*(i + (j*m))];
    }
  }
}
/*
  for (uint32 i = 0; i < A.nRows(); ++i) {
    for (uint32 j = 0; j < B.nRows(); ++j) {
#if __F4RT_DEBUG
      std::cout << i << "--" << j << std::endl;
#endif
      for (uint32 k = 0; k < B.nCols(); k++) {
#if __F4RT_DEBUG
        std::cout << "A: " << A(i,k) << " B: " << B(j,k) << std::endl;
#endif
        (*this)(i,j) += A(i,k) * B(j,k);
      }
    }
  }
*/
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

void Matrix::generateRandomMatrix(const uint32 nr, const uint32 nc, bool cmp = false) {
  m = nr;
  n = nc;
  std::ostringstream fileName;
  // generate current time for fileName
  std::string strTime;
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  strTime.append(std::to_string(now->tm_year + 1900));
  strTime.append("-");
  if (now->tm_mon+1 < 10)
    strTime.append("0");
  strTime.append(std::to_string(now->tm_mon + 1));
  strTime.append("-");
  if (now->tm_mday < 10)
    strTime.append("0");
  strTime.append(std::to_string(now->tm_mday));
  strTime.append("-");
  if (now->tm_hour < 10)
    strTime.append("0");
  strTime.append(std::to_string(now->tm_hour));
  strTime.append("-");
  if (now->tm_min < 10)
    strTime.append("0");
  strTime.append(std::to_string(now->tm_min));
  strTime.append("-");
  if (now->tm_sec < 10)
    strTime.append("0");
  strTime.append(std::to_string(now->tm_sec));

  fileName << "random-uint32-mat-" << m << "-" << n << "-" << strTime << ".mat";

  FILE* file  = fopen(fileName.str().c_str(), "wb");
  writeOne(m, file);
  writeOne(n, file);
  
  // generate a matrix of size m*n with random entries of type uint32 (but
  // being unsigned integer < 2^16 in order to get a correct multiplication)
  //std::vector<uint16> v;
  entries.resize(m*n);
  srand(time(NULL));
  for (uint64 i = 0; i < m*n; ++i) {
    entries[i] = getRandom();
  }

  writeMany(entries, file);
  fclose(file);
  if (cmp) {
    // read it from file and recheck
    Matrix B;
    file  = fopen(fileName.str().c_str(), "rb");
    B.read(file);
    fclose(file);
    check(*this, B);
    int printMatrices = 0;
    std::cout << "Matrices coincide. Print them? (1=yes, 0=no) ";
    std::cin >> printMatrices;
    if (printMatrices)
      B.print();
    B.clear();
  }

  // delete matrix
  entries.clear();
}

void Matrix::copy(const Matrix& M) {
  m = M.m;
  n = M.n;
  l = M.l;
  entries.resize(M.entries.size());
  for (uint64 i=0; i < nEntries(); ++i) 
    entries[i]  = M.entries[i];
  return;
}
