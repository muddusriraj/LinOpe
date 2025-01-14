#include <iostream>
#include <vector>
#include <cmath>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/* wrapper for std::vector<double> so that the vector can be internally used
 * in python, as python will not be able to directly access an object of type std::vector<double>*/
class darray{
    std::vector<double> content;
    public:
    darray(std::vector<double> value){
        content=value;
    }
    darray(int size){
        content.resize(size, 0);
    }
    std::vector<double> getVector() const{
        return content;
    }
    double& operator[](size_t index) {
        return content.at(index);
    }

    const double& operator[](size_t index) const {
        return content.at(index);
    }

    int getSize() const{
        return content.size();
    }

    double operator*(const darray& v) const{
        if(content.size()!=v.getSize()){
            return 0.0;
        }
        double dProd = 0;
        for(int i=0; i<content.size(); i++){
            dProd += content[i]*v[i];
        }
        return dProd;
    }

    darray operator*(double x) const{
        darray tempD(content.size());
        for(int i=0; i<content.size(); i++){
            tempD[i] = content[i]*x;
        }
        return tempD;
    }

    darray operator+(darray v) const{
        darray sum = v.getVector();
        if(content.size()!=v.getSize()){
            return darray(1);
        }
        for(int i=0; i<content.size(); i++){
            sum[i] += content[i];
        }
        return sum;
    }
    darray operator-(darray v) const{
        darray sum = v.getVector();
        if(content.size()!=v.getSize()){
            return darray(1);
        }
        for(int i=0; i<content.size(); i++){
            sum[i] -= content[i];
        }
        return sum;
    }

    double norm(){
        return sqrt((*this) * (*this));
    }

};

class Matrix{
    std::vector<darray> matrixstore; 
    int rownum;
    int colnum;
public:
    Matrix(int m, int n){
        matrixstore.resize(m, darray(std::vector<double>(n, 0.0))); 
        rownum=m;
        colnum=n;
    }
    Matrix(const Matrix& A){
        matrixstore.resize(A.rows(), darray(std::vector<double>(A.cols(), 0.0)));
        rownum=A.rows();
        colnum=A.cols();
        for(int i=0;i<rownum;i++){
            matrixstore[i] = A[i];
        }
    }
    darray& operator[](int i){
        return matrixstore[i];
    }
    const darray& operator[](int i) const{
        return matrixstore[i];
    }
    int rows() const{
        return rownum;
    }
    int cols() const{
        return colnum;
    }
    void reduce(){
        for(int i=0;i<(*this).rows();i++){
            if((*this)[i][i]==0){
                for(int j=i;j<(*this).rows();j++){
                    if((*this)[j][i]!=0){
                        darray temp = (*this)[j];
                        (*this)[j]=(*this)[i];
                        (*this)[i]=temp;
                    }
                }
            }
        }
        for(int i=0;i<(*this).rows();i++){
            for(int j=i+1;j<(*this).rows();j++){
                if(j>=(*this).rows()){
                    break;
                }
                if((*this)[j][i]!=0){
                    double factor = (*this)[j][i]/(*this)[i][i];
                    for(int x=0;x<(*this).cols();x++){
                        (*this)[j][x] = (*this)[j][x] - factor*(*this)[i][x];
                    }   
                }
            }
        }
    }
    Matrix operator* (double x){
        for(int i=0;i<colnum;i++){
            for(int j=0; j<rownum;j++){
                (*this)[i][j] = (*this)[i][j] * x;
            }
        }
        return (*this);
        
    }
    darray operator* (darray vectA){
        std::vector<double> vect = vectA.getVector();
        if(vect.size()==colnum){
            std::vector<double> solutionvect(rownum, 0);
            for(int n=0;n<colnum;n++){
                for(int m=0; m<rownum; m++){
                    solutionvect[m] = solutionvect[m] + vect[n]*(*this)[m][n];
                }
            }
            return darray(solutionvect);
        }
        throw std::invalid_argument("Invalid dimensions for Vector Matrix multiplication");
    }
    Matrix operator* (Matrix B){
        if(colnum==B.rows()){
            Matrix solutionMatrix = Matrix(rownum, B.cols());
            for(int m=0; m<rownum;m++){
                for(int x=0; x<B.cols(); x++){
                    for(int n=0; n<colnum;n++){
                        solutionMatrix[m][x] = solutionMatrix[m][x] + (*this)[m][n]*B[n][x];
                    }
                }
            }
            return solutionMatrix;
        }
        throw std::invalid_argument("Attempting to Multiply Matrices of Incorrect Dimensions");
    }


};

darray subtract(darray A, darray B){
    if(A.getSize()==B.getSize()){
        std::vector<double> vectA = A.getVector();
        std::vector<double> vectB = B.getVector();
        std::vector<double> vectC(vectA.size(), 0);
        for(int i=0;i<vectA.size();i++){
            vectC[i] = vectA[i]-vectB[i];
        }
        return darray(vectC);
    }
    throw std::invalid_argument("Invalid dimensions for vector subtraction");
}

darray add(darray A, darray B){
    if(A.getSize()==B.getSize()){
        std::vector<double> vectA = A.getVector();
        std::vector<double> vectB = B.getVector();
        std::vector<double> vectC(vectA.size(), 0);
        for(int i=0;i<vectA.size();i++){
            vectC[i] = vectA[i]+vectB[i];
        }
        return darray(vectC);
    }
    throw std::invalid_argument("Invalid dimensions for vector addition");
}

Matrix add(Matrix A, Matrix B){
    if(A.rows()==B.rows() && A.cols()==B.cols()) {
        Matrix SolutionMatrix(A.rows(), A.cols());
        for(int m=0;m<A.rows();m++){
            for(int n=0;n<A.cols();n++){
                SolutionMatrix[m][n]=A[m][n]+B[m][n];
            }
        }
        return SolutionMatrix;
    }
    throw std::invalid_argument("Attempting to Add Matrices of Different Dimensions");
}

Matrix subtract(Matrix A, Matrix B){
    if(A.rows()==B.rows() && A.cols()==B.cols()) {
        Matrix SolutionMatrix(A.rows(), A.cols());
        for(int m=0;m<A.rows();m++){
            for(int n=0;n<A.cols();n++){
                SolutionMatrix[m][n]=A[m][n]-B[m][n];
            }
        }
        return SolutionMatrix;
    }
    throw std::invalid_argument("Attempting to Add Matrices of Different Dimensions");
}

darray multiply(double scalar, darray vectA){
    darray temp = darray(vectA.getVector());
    for(int i=0; i<vectA.getSize();i++){
        temp[i] *= scalar;
    }
    return temp;
}

Matrix multiply(double scalar, Matrix A){
    Matrix SolutionMatrix(A.rows(), A.cols());
    for(int m=0;m<A.rows();m++){
        for(int n=0;n<A.cols();n++){
            SolutionMatrix[m][n]=scalar*A[m][n];
        }
    }
    return SolutionMatrix;
};

Matrix multiply(Matrix A, Matrix B){
    if(A.cols()==B.rows()){
        Matrix solutionMatrix = Matrix(A.rows(), B.cols());
        for(int m=0; m<A.rows();m++){
            for(int x=0; x<B.cols(); x++){
                for(int n=0; n<A.cols();n++){
                    solutionMatrix[m][x] = solutionMatrix[m][x] + A[m][n]*B[n][x];
                }
            }
        }
        return solutionMatrix;
    }
    throw std::invalid_argument("Attempting to Multiply Matrices of Incorrect Dimensions");
}

darray multiply(Matrix A, darray vectA){
    std::vector<double> vect = vectA.getVector();
    if(vect.size()==A.cols()){
        std::vector<double> solutionvect(A.rows(), 0);
        for(int n=0;n<A.cols();n++){
            for(int m=0; m<A.rows(); m++){
                solutionvect[m] = solutionvect[m] + vect[n]*A[m][n];
            }
        }
        return darray(solutionvect);
    }
    throw std::invalid_argument("Invalid dimensions for Vector Matrix multiplication");
}

Matrix reduce(Matrix A){
    for(int i=0;i<A.rows();i++){
        if(A[i][i]==0){
            for(int j=i;j<A.rows();j++){
                if(A[j][i]!=0){
                    darray temp = A[j];
                    A[j]=A[i];
                    A[i]=temp;
                }
            }
        }
    }
    for(int i=0;i<A.rows();i++){
        for(int j=i+1;j<A.rows();j++){
            if(j>=A.rows()){
                break;
            }
            if(A[j][i]!=0){
                double factor = A[j][i]/A[i][i];
                for(int x=0;x<A.cols();x++){
                    A[j][x] = A[j][x] - factor*A[i][x];
                }   
            }
        }
    }
    return A;
}

Matrix rowReduce(Matrix A){
    for(int i=0;i<A.rows();i++){
        if(A[i][i]==0){
            for(int j=i;j<A.rows();j++){
                if(A[j][i]!=0){
                    darray temp = A[j];
                    A[j]=A[i];
                    A[i]=temp;
                }
            }
        }
    }
    for(int i=0;i<A.rows(); i++){
        for(int j=i+1;j<A.rows();j++){
            if(j>=A.rows()){
                break;
            }
            if(A[j][i]!=0){
                double factor = A[j][i]/A[i][i];
                for(int x=0;x<A.cols();x++){
                    A[j][x] = A[j][x] - factor*A[i][x];
                }   
            }
        }
    }
    double factor{};
    for(int i=A.rows()-1;i>=0;i--){
        if(A[i][i]==0) continue;
        std::cout<<factor;
        factor = 1.0/A[i][i];
        A[i] = A[i] * factor;


    }
    return A;
}

Matrix transpose(Matrix A){
    Matrix transposeA(A.cols(), A.rows());
    for(int i=0; i<A.cols(); i++){
        for(int j=0; j<A.rows(); j++){
            transposeA[j][i]=A[i][j];
        }
    }
    return transposeA;
}

Matrix identity(int x){
    Matrix A(x, x);
    for(int i=0;i<x;i++){
        A[i][i]=1;
    }
    return A;
}

Matrix inverse(Matrix A){
    Matrix invMatr = identity(A.cols());
    return invMatr;

}

//determinant
double det(Matrix A){
    double coeff=1;
    double det=1;
    for(int i=0;i<A.rows();i++){
        if(A[i][i]==0){
            for(int j=i;j<A.rows();j++){
                if(A[j][i]!=0){
                    coeff*=-1;
                    darray temp = A[j];
                    A[j]=A[i];
                    A[i]=temp;
                }
            }
        }
    }
    for(int i=0;i<A.rows();i++){
        for(int j=i+1;j<A.rows();j++){
            if(j>=A.rows()){
                break;
            }
            if(A[j][i]!=0){
                double factor = A[j][i]/A[i][i];
                for(int x=0;x<A.cols();x++){
                    A[j][x] = A[j][x] - factor*A[i][x];
                }   
            }
        }
    }
    if(A.rows()==A.cols()){
        for(int i=0;i<A.rows(); i++){
            det *= A[i][i];
        }
        det*=coeff;

    }
    return det;
}

darray normalize(darray A){
    darray vectAns = darray(A);
    double ans = sqrt(A*A);
    for(int i=0; i<A.getSize();i++){
        vectAns[i] = vectAns[i]/ans;
    }
    return vectAns;
}

void printMatrix(Matrix A){
    for(int i=0; i<A.rows(); i++){
        for(int j=0;j<A.cols();j++){
            std::cout << A[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void printVector(darray A){
    std::vector<double> temp = A.getVector();
    for(int i=0;i<temp.size();i++){
        std::cout << temp[i] << " ";
    }
}

PYBIND11_MODULE(operators, m) {
    py::class_<darray>(m, "darray")
        .def(py::init<std::vector<double>>())
        .def(py::init<int>())
        .def("getVector", &darray::getVector)
        //Had to create lambdas for certain methods due to errors when attempting to overload cast
        .def("__getitem__", [](const darray& self, size_t index) -> double {
            return self[index];
        }, py::arg("index"))
        .def("__setitem__", [](darray& self, size_t index, double value){
            self[index]=value;
        }, py::arg("index"), py::arg("value"))
        .def("getSize", &darray::getSize)
        .def("__mul__", [](const darray& self, darray v) -> double {
            return self*v;
        }, py::arg("v"))
        .def("__mul__", [](const darray&self, double x) -> darray{
            return self*x;
        }, py::arg("x"))
        .def("__add__", &darray::operator+, py::arg("v"))
        .def("__sub__", &darray::operator-, py::arg("v"))
        .def("norm", &darray::norm);

    m.def("normalize", &normalize);  

    py::class_<Matrix>(m, "Matrix")
        .def(py::init<int, int>())
        .def("__getitem__", [](const Matrix& self, int i){
            return self[i];
        }, py::arg("i"))
        .def("__setitem__", [](Matrix& self, int i, darray value){
            self[i] = value;
        }, py::arg("i"), py::arg("value"))
        .def("rows", &Matrix::rows)
        .def("cols", &Matrix::cols)
        .def("reduce", &Matrix::reduce)
        .def("__mul__", py::overload_cast<double>(&Matrix::operator*), py::arg("x"))
        .def("__mul__", py::overload_cast<darray>(&Matrix::operator*), py::arg("vectA"))
        .def("__mul__", py::overload_cast<Matrix>(&Matrix::operator*), py::arg("B"));

    m.def("subtract", py::overload_cast<darray, darray>(&subtract), "Subtract two vectors");
    m.def("add", py::overload_cast<darray, darray>(&add), "Add two vectors");
    m.def("add", py::overload_cast<Matrix, Matrix>(&add), "Add two matrices");
    m.def("subtract", py::overload_cast<Matrix, Matrix>(&subtract), "Subtract two matrices");
    m.def("multiply", py::overload_cast<double, darray>(&multiply), "Multiply a vector by a scalar");
    m.def("multiply", py::overload_cast<double, Matrix>(&multiply), "Multiply a matrix by a scalar");
    m.def("multiply", py::overload_cast<Matrix, Matrix>(&multiply), "Multiply two matrices");
    m.def("multiply", py::overload_cast<Matrix, darray>(&multiply), "Multiply a matrix by a vector");
    m.def("reduce", py::overload_cast<Matrix>(&reduce), "Reduce a matrix to row echelon form");
    m.def("row_reduce", py::overload_cast<Matrix>(&rowReduce), "Row reduce a matrix");
    m.def("transpose", py::overload_cast<Matrix>(&transpose), "Transpose a matrix");
    m.def("identity", py::overload_cast<int>(&identity), "Generate an identity matrix of size x");
    m.def("inverse", py::overload_cast<Matrix>(&inverse), "Find the inverse of a matrix");
    m.def("det", py::overload_cast<Matrix>(&det), "Calculate the determinant of a matrix");
    m.def("normalize", py::overload_cast<darray>(&normalize), "Normalize a vector");
    m.def("print_matrix", py::overload_cast<Matrix>(&printMatrix), "Print a matrix");
    m.def("print_vector", py::overload_cast<darray>(&printVector), "Print a vector");

}
