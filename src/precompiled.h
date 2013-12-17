#ifndef PRECOMPILED_H
#define PRECOMPILED_H

#if defined __cplusplus

#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <string>
#include <iomanip>
#include <cmath>
#include <assert.h>
#include <cstdlib>
#include <pmmintrin.h>
#include <memory>
#include <cstring>

#ifdef WITH_OMP
#include <omp.h>
#endif

typedef unsigned long long ull;
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;

using std::shared_ptr;
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::istreambuf_iterator;

#endif

#endif // PRECOMPILED_H
