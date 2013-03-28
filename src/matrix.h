/**
 * \file   matrix.h
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Header file for matrix implementation.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

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
#include <typeinfo>
#include <sys/time.h>
#include <cassert>
#ifdef __F4RT_HAVE_INTEL_TBB
#include <tbb/tbb.h>
#endif
#ifdef __F4RT_HAVE_PTHREAD_H
#include <pthread.h>
#endif
#ifdef __F4RT_HAVE_OPENMP
#include <omp.h>
#endif
#if defined(__F4RT_HAVE_KAAPI)
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


// matrix entry type:
// if uint16 entries are used we store uint64
// if float entries are used we store double
typedef uint64 mat; // mat entry type
typedef uint16 rmat; // real type of entry

// Dense matrix
struct matrix {
  uint32 rows;
  uint32 cols;
  mat *entries;
};



void clear(matrix& M);
void read(matrix& M, FILE* file);
void write(const matrix M, FILE* file);
void copy(const matrix& M, matrix& N);
void print(const matrix& M);
void transpose(const matrix& M, matrix& N);
rmat getRandVal();
void genRandom( const uint32 m, const uint32 n, 
                bool cmp, bool timestamp);



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
}
#endif
