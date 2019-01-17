/*=========================================================================
** Title       : Interface for matrixarith.c
** File        : matrixarith.h
** =======================================================================*/

#ifndef __MATRIXARITH__
#define __MATRIXARITH__

#include "mexheader.h"

void matrixMultiply(mxArrayConstPtr R, mxArrayConstPtr A, mxArrayConstPtr B);
void matrixSubtract(mxArrayConstPtr R, mxArrayConstPtr A, mxArrayConstPtr B);
double matrixFrobeniusNormSquared(mxArrayConstPtr A);

#endif
