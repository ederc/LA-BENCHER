#include <map>
#include <cmath>
#include <algorithm>
#include <set>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <time.h>
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
  uint32 m;
  uint32 n;
  uint64 l;
  std::vector<uint16> entries;
  
  void read(FILE* file);
  void write(FILE* file);
public:
  Matrix() : m(0), n(0) {}
  Matrix(size_t m_, size_t n_) : m(m_), n(n_), l(m_*n_), entries(m_*n_) {}
  Matrix(size_t m_, size_t n_, uint16 val) : m(m_), n(n_), l(m_*n_), entries(m_*n_, val) {}

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

  uint16* column(size_t j) {
    return &*(entries.begin() + j*m);
  } 
  
  uint16& operator()(size_t i, size_t j)  {
    return entries[i + (j * m)];
  }
  
  const uint16& operator()(size_t i, size_t j) const  {
    return entries[i + (j * m)];
  }

  void generateRandomMatrix(const uint32 m, const uint32 n, bool cmp);

  void print();
};
