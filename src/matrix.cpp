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

/*
// matrix deletion
void clear(matrix& M) {
  M.rows  = 0;
  M.cols  = 0;
  M.entries.clear();
  M.entries = NULL;
}

// matrix reading
// memory for M is allocated inside!
void read(matrix& M, FILE* file) {
  if (fread(&M.rows, sizeof(uint32), 1, file) != 1)
    std::cerr << "Error while reading file." << std::endl;
  if (fread(&M.cols, sizeof(uint32), 1, file) != 1)
    std::cerr << "Error while reading file." << std::endl;

  M.entries = (mat *)malloc(M.rows * M.cols * sizeof(mat));
  if (fread(M.entries, sizeof(mat), sizeof(M.entries)/sizeof(mat), file)
      != sizeof(M.entries)/sizeof(mat))
    std::cerr << "Error while reading file." << std::endl;
}

// matrix writing
void write(const matrix M, FILE* file) {
  if (fwrite(&M.rows, sizeof(uint32), 1, file) != 1)
    std::cerr << "Error while writing to file." << std::endl;
  if (fwrite(&M.cols, sizeof(uint32), 1, file) != 1)
    std::cerr << "Error while writing to file." << std::endl;
  if (M.entries == NULL)
    return;
  if (fwrite(M.entries, sizeof(mat), sizeof(M.entries)/sizeof(mat), file) 
      != sizeof(M.entries)/sizeof(mat))
    std::cerr << "Error while writing to file." << std::endl;
}

// matrix copying
void copy(const matrix& M, matrix& N) {
  N.rows    = M.rows;
  N.cols    = M.cols;
  N.entries = (mat *)malloc(N.rows * N.cols * sizeof(mat));
  for (uint32 i = 0; i < N.rows * N.cols; i++)
    N.entries[i]  = M.entries[i];
}

// matrix printing
void print(const matrix& M) {
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
void transpose(const matrix& M, matrix& N) {
  uint32 tempRows     = M.cols;
  uint32 tempCols     = M.rows;
  mat *tempEntries = (mat *)malloc(tempRows * tempCols * sizeof(mat));
  uint32 i,j;
  for (i = 0; i < M.rows; ++i) {
    for (j = 0; j < M.cols; ++j) {
      tempEntries[i+j*M.rows]  = M.entries[j+i*M.cols];
    }
  }
  N.rows    = tempRows;
  N.cols    = tempCols;
  N.entries = tempEntries;
}

// get random entries
rmat getRandVal() {
  return static_cast<rmat>(std::rand());
}

// random matrix generating
void genRandom( const uint32 m, const uint32 n, 
                bool cmp, bool timestamp) {
  struct matrix M;
  M.rows     = m;
  M.cols     = n;
  std::cout << m << " -- " << n << std::endl;
  std::cout << M.rows << " -- " << M.cols << std::endl;
  std::cout << sizeof(M.entries) << std::endl;
  std::cout << sizeof(mat) << std::endl;
  std::cout << M.rows * M.cols * sizeof(*M.entries) << std::endl;
  mat *temp   = (mat *)malloc(M.rows * M.cols * sizeof(*M.entries));
  M.entries   = temp; 
  std::cout << sizeof(M.entries) << std::endl;
  std::cout << sizeof(M.entries[0]) << std::endl;
 
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

  fileName << "random-mat-" << m << "-" << n << "-" << strTime << ".mat";

  FILE* file  = fopen(fileName.str().c_str(), "ab");
  srand(time(NULL));
  for (uint64 i = 0; i < M.rows * M.cols; ++i) {
    M.entries[i] = getRandVal();
    std::cout << M.entries[i] << std::endl;
  }
  std::cout << sizeof(temp) << std::endl;
  std::cout << sizeof(temp[0]) << std::endl;
  write(M, file); 
  fclose(file);

  if (cmp) {
    // read it from file and recheck
    matrix N;
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

*/
////////////////////////

// not coincide: return 1
// coincide:     return 0
int check(const Matrix& A, const Matrix& B, int unittest = 0) {
  if (!unittest)
    std::cout << "Checking generated matrix with data written to file." << std::endl;
  for (uint32 i = 0; i < A.nRows(); ++i)
    for (uint32 j = 0; j < A.nCols(); ++j) 
      if (A(i,j) != B(i,j)) {
        if (!unittest)
          std::cerr << "Error: Matrices do not coincide at position (" << i << "," << j << ")" << std::endl;
        return 1;
      }
  return 0;
}
void Matrix::read(FILE* file) {
  const uint32 rowCount   = readOne(file);
  const uint32 colCount   = readOne(file);
  m = rowCount;
  n = colCount;
  uint64 dim  = static_cast<uint64>(rowCount)*colCount;
  entries.resize(dim);
  readMany(file, dim, entries);
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

mat getRandomVal() {
  return static_cast<rmat>(std::rand());
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
    char buffer[50];
    sprintf(buffer,"%d",now->tm_year + 1900);
    strTime.append(buffer);
    strTime.append("-");
    if (now->tm_mon+1 < 10)
      strTime.append("0");
    sprintf(buffer,"%d",now->tm_mon + 1);
    strTime.append(buffer);
    strTime.append("-");
    if (now->tm_mday < 10)
      strTime.append("0");
    sprintf(buffer,"%d",now->tm_mday);
    strTime.append(buffer);
    strTime.append("-");
    if (now->tm_hour < 10)
      strTime.append("0");
    sprintf(buffer,"%d",now->tm_hour);
    strTime.append(buffer);
    strTime.append("-");
    if (now->tm_min < 10)
      strTime.append("0");
    sprintf(buffer,"%d",now->tm_min);
    strTime.append(buffer);
    strTime.append("-");
    if (now->tm_sec < 10)
      strTime.append("0");
    sprintf(buffer,"%d",now->tm_sec);
    strTime.append(buffer);

    fileName << "random-mat-" << m << "-" << n << "-" << strTime << ".mat";
  } else {
    fileName << "random-mat-" << m << "-" << n << ".mat";
  }

  FILE* file  = fopen(fileName.str().c_str(), "ab");
  writeOne(m, file);
  writeOne(n, file);
  
  // generate a matrix of size m*n with random entries of type uint32 (but
  // being unsigned integer < 2^16 in order to get a correct multiplication)
  //std::vector<uint16> v;
  uint64 dim  = static_cast<uint64>(m)*n;
  entries.resize(dim);
  srand(time(NULL));
  for (uint64 i = 0; i < dim; ++i) {
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
  entries.resize(M.entries.size());
  for (uint64 i=0; i < m * n; ++i) 
    entries[i]  = M.entries[i];
  return;
}

// self-transpose
void Matrix::transpose() {
  uint32 tempRows = n;
  uint32 tempCols = m;
  std::vector<mat> tempEntries;
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
  entries.resize(M.entries.size());
  for (uint32 i=0; i<M.nRows(); ++i) {
    for (uint32 j=0; j<M.nCols(); ++j) {
      (*this)(j,i)  = M(i,j);
    }
  }
  return;
}
