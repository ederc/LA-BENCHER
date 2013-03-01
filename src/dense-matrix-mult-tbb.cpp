#include <vector>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <iterator>
#include <cstdio>
#include <cassert>
#include <tbb/mutex.h>

enum ArgumentType {
    TYPE_NONE, TYPE_INT, TYPE_INTEGER, TYPE_DOUBLE, TYPE_STRING
};
#define TYPE_BOOL TYPE_NONE

struct Argument {
  char             c;
  const char      *example;
  const char      *helpString;
  ArgumentType     type;
  void            *data;
};

class DenseRow {
public:
  typedef uint16_t Scalar;
  typedef uint32_t ScalarProduct;
  typedef uint64_t ScalarProductSum;

  static ScalarProduct product(const Scalar a, const Scalar b) {
    return static_cast<ScalarProduct>(a) * b;
  }

  static void multiplyAdd(
    const Scalar a,
    const Scalar b,
    ScalarProductSum& x
  ) {
    x += product(a, b);
  }

  static void add(const Scalar a, ScalarProductSum& sum) {
    std::cout << "ADDING: "<< sum << " + " << a << std::endl;
    sum += a;
  }

  Scalar modulusOf(ScalarProductSum x, Scalar modulus) {
    return static_cast<Scalar>(x % modulus);
  }

  DenseRow() {}
  DenseRow(size_t colCount): mEntries(colCount) {}


  size_t colCount() const {return mEntries.size();}
  bool empty() const {return mEntries.empty();}

  void clear(size_t colCount = 0) {
    mEntries.clear();
    mEntries.resize(colCount);
  }

  ScalarProductSum& operator[](size_t col) {
    assert(col < colCount());
    return mEntries[col];
  }

  ScalarProductSum const& operator[](size_t col) const {
    assert(col < colCount());
    return mEntries[col];
  }

private:
  std::vector<ScalarProductSum> mEntries;
};

int main(int argc, char **argv)
{
	const char *fileName = "";

	bool pass = true;
	bool validate_results = false;
	bool free_mem = false;
	bool compute_Rref = false;
	bool use_standard_method = false;
	int n_threads = 8;
	bool horizontal = false;
	bool reconstruct_old = false;

	static Argument args[] =
	{
		{ 'f', "-f File", "The file name where the matrix is stored", TYPE_STRING, &fileName },
		{ 'p', "-p NUM_THREADS", "Number of threads (DEFAULT 8)", TYPE_INT, &n_threads },
		{ '\0' }
	};

	parseArguments(argc, argv, args, "", 0);
	free_mem = !free_mem;
  std::cout << "test" << std::endl;

}

