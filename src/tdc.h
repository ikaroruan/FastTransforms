#ifndef TDC_H
#define TDC_H

/*
For assignments of the form:
for (int i = i.start; i < i.stop; i++)
    x[i] = ;
*/
typedef struct {
    int start;
    int stop;
} unitrange;

#define CONCAT(prefix, name, suffix) prefix ## name ## suffix

#define hash(m,n) hash[(m)+(n)*M]
#define hierarchicalmatrices(m,n) hierarchicalmatrices[(m)+(n)*M]
#define densematrices(m,n) densematrices[(m)+(n)*M]
#define lowrankmatrices(m,n) lowrankmatrices[(m)+(n)*M]
#define fmmblocks(m,n) fmmblocks[(m)+(n)*M]

#define FLT float
#define X(name) CONCAT(, name, f)
#include "tridiagonal_source.h"
#include "triangular_banded_source.h"
#include "hierarchical_source.h"
#include "dprk_source.h"
#include "tdc_source.h"
#undef FLT
#undef X

#define FLT double
#define X(name) CONCAT(, name, )
#include "tridiagonal_source.h"
#include "triangular_banded_source.h"
#include "hierarchical_source.h"
#include "dprk_source.h"
#include "tdc_source.h"
#undef FLT
#undef X

#define FLT long double
#define X(name) CONCAT(, name, l)
#include "tridiagonal_source.h"
#include "triangular_banded_source.h"
#include "hierarchical_source.h"
#include "dprk_source.h"
#include "tdc_source.h"
#undef FLT
#undef X

#define FLT __float128
#define X(name) CONCAT(, name, q)
#include "tridiagonal_source.h"
#include "triangular_banded_source.h"
#include "hierarchical_source.h"
#include "dprk_source.h"
#include "tdc_source.h"
#undef FLT
#undef X

#endif // TDC_H