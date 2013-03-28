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

// Dense floating point entry matrix
template <class T>
struct matrix {
  uint32 rows;
  uint32 cols;
  T *entries;
};

// matrix deletion
template <class T>
void clear(matrix<T>& M) {
  M.rows  = 0;
  M.cols  = 0;
  free(M.entries);
  M.entries = NULL;
}

// matrix reading
// memory for M is allocated inside!
template <class T>
void read(matrix<T>& M, FILE* file) {
  if (fread(&M.rows, sizeof(uint32), 1, file) != 1)
    std::cerr << "Error while reading file." << std::endl;
  if (fread(&M.cols, sizeof(uint32), 1, file) != 1)
    std::cerr << "Error while reading file." << std::endl;

  M.entries = (T *)malloc(M.rows * M.cols * sizeof(T));
  if (fread(M.entries, sizeof(T), sizeof(M.entries)/sizeof(T), file)
      != sizeof(M.entries)/sizeof(T))
    std::cerr << "Error while reading file." << std::endl;
}

// matrix writing
template <class T>
void write(const matrix<T> M, FILE* file) {
  if (fwrite(&M.rows, sizeof(uint32), 1, file) != 1)
    std::cerr << "Error while writing to file." << std::endl;
  if (fwrite(&M.cols, sizeof(uint32), 1, file) != 1)
    std::cerr << "Error while writing to file." << std::endl;
  if (M.entries == NULL)
    return;
  if (fwrite(M.entries, sizeof(T), sizeof(M.entries)/sizeof(T), file) 
      != sizeof(M.entries)/sizeof(T))
    std::cerr << "Error while writing to file." << std::endl;
}

// matrix copying
template <class T>
matrix<T> *copy(const matrix<T>& M) {
  matrix<T> N;
  N.rows    = M.rows;
  N.cols    = M.cols;
  N.entries = (T *)malloc(N.rows * N.cols * sizeof(T));
  for (int i = 0; i < N.rows * N.cols; i++)
    N.entries[i]  = M.entries[i];

  return N;
}

// matrix printing
template <class T>
void print(const matrix<T>& M) {
  std::cout << "#rows" << std::setw(15) << M.rows << std::endl;
  std::cout << "#cols" << std::setw(15) << M.cols << std::endl;
  if (M.rows * M.cols > 1000) {
    int pr;
    std::cout << "NOTE: The matrix consists of " << M.rows * M.cols << " entries." << std::endl;
    std::cout << "      Do you really want to print it? (1=yes, else=no) ";
    std::cin >> pr;
    if (pr != 1)
      return;
  }
  uint32 i,j;
  for (i = 0; i < M.rows; ++i) {
    for (j = 0; j < M.cols-1; ++j) {
      std::cout  << std::setw(5) << M.entries[j+i*M.cols] << " | ";
    }
    std::cout  << std::setw(5) << M.entries[j+i*M.cols];
    std::cout << std::endl;
  }
}

// matrix transposing
template <class T>
matrix<T> *transpose(const matrix<T>& M) {
  matrix<T> N;
  uint32 tempRows     = M.cols;
  uint32 tempCols     = M.rows;
  uint32 *tempEntries = (T *)malloc(tempRows * tempCols * sizeof(T));
  uint32 i,j;
  for (i = 0; i < M.rows; ++i) {
    for (j = 0; j < M.cols; ++j) {
      tempEntries[i+j*M.rows]  = M.entries[j+i*M.cols];
    }
  }
  N.rows    = tempRows;
  N.cols    = tempCols;
  N.entries = tempEntries;

  return &N;
}

// get random entries
// T t only used for template type deduction
template <typename T>
T getRandVal(T t) {
  return static_cast<T>(t * std::rand());
}


// random matrix generating
template <class T>
void genRandom( matrix<T> M, const uint32 m, const uint32 n, 
                const uint32 type, bool cmp, bool timestamp) {
  M.rows     = m;
  M.cols     = n;
  M.entries  = (T *)malloc(M.rows * M.cols * sizeof(T));
 
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

  // uint16
  if (type == 1) {
    fileName << "random-uint16-mat-" << m << "-" << n << "-" << strTime << ".mat";
  }
  // uint32
  if (type == 2) {
    fileName << "random-uint32-mat-" << m << "-" << n << "-" << strTime << ".mat";
  }
  // float
  if (type == 3) {
    fileName << "random-float-mat-" << m << "-" << n << "-" << strTime << ".mat";
  }
  FILE* file  = fopen(fileName.str().c_str(), "ab");
  srand(time(NULL));
  T temp = 1;
  for (uint64 i = 0; i < M.rows * M.cols; ++i) {
    M.entries[i] = getRandVal(temp);
  }

  write(M, file); 
  fclose(file);

  if (cmp) {
    // read it from file and recheck
    matrix<T> N;
    file  = fopen(fileName.str().c_str(), "rb");
    read(N, file);
    fclose(file);
    //TOOOOOOO DOOOOOO
    //check(*this, B);
    int printMatrices = 0;
    std::cout << "Matrices coincide. Print them? (1=yes, 0=no) ";
    std::cin >> printMatrices;
    if (printMatrices)
      print(N);
    clear(N);
  }

  clear(M);
}

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
