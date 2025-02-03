#ifndef LRM_H
#define LRM_H

#include "graph.h"

darray& lrm(Matrix &X, Matrix &Y);

std::vector<double> slrm(darray x, darray y);

#endif