#include "lrm.h"

darray& lrm(Matrix &X, Matrix &Y) {
	darray beta = transpose(inverse(transpose(X) * X) * (transpose(X) * Y))[0];

	return beta;

}

std::vector<double> slrm(darray x, darray y) {
	double beta = (x * y) / (x * x);
	double intercept = (y.sum() / y.getSize()) - (beta * (x.sum() / x.getSize()));

	return { beta, intercept };

}