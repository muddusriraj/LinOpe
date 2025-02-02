#ifndef OPERATORS_H
#define OPERATORS_H

#include <iostream>
#include <vector>
#include <cmath>

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>

// namespace py = pybind11;

class darray{
    std::vector<double> content;
    public:
    darray(std::vector<double> value);

    darray(int size);

    std::vector<double> getVector() const;
    
    double& operator[](size_t index);

    // Const operator[] for getting values
    const double& operator[](size_t index) const;

    int getSize() const;

    double operator*(const darray v) const;

    darray operator*(double x) const;

    darray operator+(darray v) const;

    darray operator-(darray v) const;

    double norm();

    double sum();

};

class Matrix{
    std::vector<darray> matrixstore; 
    int rownum;
    int colnum;
public:
    Matrix(int m, int n);

    Matrix(const Matrix& A);

    darray& operator[](int i);

    const darray& operator[](int i) const;
    
    int rows() const;

    int cols() const;

    void reduce();

    Matrix operator* (double x);

    darray operator* (darray vectA);

    Matrix operator* (Matrix A);

};

darray subtract(darray A, darray B);

darray add(darray A, darray B);

Matrix add(Matrix A, Matrix B);

Matrix subtract(Matrix A, Matrix B);

darray multiply(double scalar, darray vectA);

Matrix multiply(double scalar, Matrix A);

Matrix multiply(Matrix A, Matrix B);

darray multiply(Matrix A, darray vectA);

Matrix reduce(Matrix A);

Matrix rowReduce(Matrix A);

Matrix transpose(Matrix A);

Matrix identity(int x);

double det(Matrix A);

darray normalize(darray A);

void printMatrix(Matrix A);

void printVector(darray A);

Matrix inverse(Matrix A);

double mean(Matrix A);

#endif