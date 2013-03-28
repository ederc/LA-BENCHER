/**
 * \file   matrix.cpp
 * \author Christian Eder ( christian.eder@inria.fr )
 * \date   March 2013
 * \brief  Source file for matrix implementation.
 *         This file is part of F4RT, licensed under the GNU General
 *         Public License version 3. See COPYING for more information.
 */

#include "matrix.h"

#define __F4RT_DEBUG  0
////////////////////////

// not coincide: return 1
// coincide:     return 0
int check(const Matrix& A, const Matrix& B, int unittest = 0) {
  if (!unittest)
    std::cout << "Checking generated matrix with data written to file." << std::endl;
  for (size_t i = 0; i < A.nRows(); ++i)
    for (size_t j = 0; j < A.nCols(); ++j) 
      if (A(i,j) != B(i,j)) {
        if (!unittest)
          std::cerr << "Error: Matrices do not coincide at position (" << i << "," << j << ")" << std::endl;
        return 1;
      }
  return 0;
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
  std::cout << "#rows" << std::setw(15) << m << std::endl;
  std::cout << "#cols" << std::setw(15) << n << std::endl;
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
      std::cout  << std::setw(5) << entries[j+i*n] << " | ";
    }
    std::cout  << std::setw(5) << entries[j+i*n];
    std::cout << std::endl;
  }
}

float getRandomVal() {
  return static_cast<float>(std::rand());
}


void Matrix::generateRandomMatrix(const uint32 nr, const uint32 nc, bool cmp = false, 
                                  bool timestamp = false) {
  m = nr;
  n = nc;
  std::ostringstream fileName;
  // generate current time for fileName
  if (timestamp) {
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

    fileName << "random-float-mat-" << m << "-" << n << "-" << strTime << ".mat";
  } else {
    fileName << "random-float-mat-" << m << "-" << n << ".mat";
  }

  FILE* file  = fopen(fileName.str().c_str(), "ab");
  writeOne(m, file);
  writeOne(n, file);
  
  // generate a matrix of size m*n with random entries of type uint32 (but
  // being unsigned integer < 2^16 in order to get a correct multiplication)
  //std::vector<uint16> v;
  entries.resize(m*n);
  srand(time(NULL));
  for (uint64 i = 0; i < m*n; ++i) {
    entries[i] = getRandomVal();
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

// self-transpose
void Matrix::transpose() {
  uint32 tempRows = n;
  uint32 tempCols = m;
  std::vector<float> tempEntries;
  tempEntries.resize(entries.size());
  for (uint32 i=0; i<m; ++i) {
    for (uint32 j=0; j<n; ++j) {
      //tempEntries[i+j*m]  = (*this)(i,j);
      tempEntries[i+j*m]  = entries[j+i*n];
    }
  }
  m = tempRows;
  n = tempCols;
  entries.clear();
  entries  = tempEntries;
  return;
}
void Matrix::transpose(const Matrix& M) {
  m = M.n;
  n = M.m;
  l = M.l;
  entries.resize(M.entries.size());
  for (uint32 i=0; i<M.nRows(); ++i) {
    for (uint32 j=0; j<M.nCols(); ++j) {
      (*this)(j,i)  = M(i,j);
    }
  }
  return;
}
