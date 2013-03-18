#ifndef F4RT_MATRIX_H
#define F4RT_MATRIX_H

#include <f4rt-config.h>
#include <map>
#include <cmath>
#include <algorithm>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <cassert>
#ifdef __F4RT_HAVE_INTEL_TBB
#include <tbb/tbb.h>
#endif
#include <pthread.h>
#ifdef __F4RT_HAVE_OPENMP
#include <omp.h>
#endif
#if defined(__F4RT_HAVE_KAAPI) && defined(__F4RT_ENABLE_KAAPI)
#include <kaapi.h>
#include <kaapic.h>
#include <kaapi++>
#endif
//#include "auxiliary.h"

typedef unsigned long long uint64;
typedef unsigned int uint32;
typedef unsigned short uint16;
typedef unsigned char uint8;

typedef signed long long int64;
typedef signed int int32;
typedef signed short int16;
typedef signed char int8;

/// Dense matrix
class Matrix {
public:
  uint32 m;
  uint32 n;
  uint64 l;
  std::vector<float> entries;
  
  void write(FILE* file);
  Matrix() : m(0), n(0) {}
  Matrix(size_t m_, size_t n_) : m(m_), n(n_), l(m_*n_), entries(m_*n_) {}
  Matrix(size_t m_, size_t n_, float val) : m(m_), n(n_), l(m_*n_), entries(m_*n_, val) {}

  void clear() {
    m = n = l = 0;
    entries.clear();
  }

  size_t nRows() const {
    return m;
  }
  
  size_t nCols() const {
    return n;
  }

  size_t nEntries() const {
    return l;
  }

  float* column(size_t j) {
    return &*(entries.begin() + j*m);
  } 
  
  float& operator()(size_t i, size_t j)  {
    return entries[j+ (i*n)];
  }
  
  const float& operator()(size_t i, size_t j) const  {
    return entries[j + (i * n)];
  }

  void resize(uint64 val) {
    entries.resize(val);
  }

  void generateRandomMatrix(const uint32 m, const uint32 n, bool cmp, bool timestamp);

  void copy(const Matrix& M);

  void transpose();

  void transpose(const Matrix& M);

  void read(FILE* file);

  void print();
};

int check(const Matrix& A, const Matrix& B, int unittest);

float getRandom();
#endif
