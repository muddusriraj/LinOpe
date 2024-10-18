#include <iostream>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
//Vector class
class darray{
    std::vector<double> content;
    public:
    darray(std::vector<double> value){
        content=value;
    }
    const std::vector<double>& getVector(){
        return content;
    }
    double& operator[](int index) {
        return content[index]; // Returns a reference to the element
    }

    // Const operator[] for getting values
    const double& operator[](int index) const {
        return content[index];
    }

    int getSize(){
        return content.size();
    }

};

//Matrix class
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
    darray& operator[](int i){
        return matrixstore[i];
    }
    const darray& operator[](int i) const{
        return matrixstore[i];
    }
    int rows(){
        return rownum;
    }
    int cols(){
        return colnum;
    }
};

//darray (C++ vector) subtraction
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

//darray (C++ vector) addition
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

//Matrix addition
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

//Matrix subtraction
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

//Matrix multiplication for scalars
Matrix multiply(double scalar, Matrix A){
    Matrix SolutionMatrix(A.rows(), A.cols());
    for(int m=0;m<A.rows();m++){
        for(int n=0;n<A.cols();n++){
            SolutionMatrix[m][n]=scalar*A[m][n];
        }
    }
    return SolutionMatrix;
};

//Multiplies matrices
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

//Matrix multiplication for vectors
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


//Solve Matrix
Matrix reduce(Matrix A){
    bool cont = true;
    std::vector<double> zeroV;
    for(int i=0;i<A.cols();i++){
        zeroV.push_back(0);
    }
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

//Prints Matrix in a readable format
void printMatrix(Matrix A){
    for(int i=0; i<A.rows(); i++){
        for(int j=0;j<A.cols();j++){
            std::cout << A[i][j] << " ";
        }
        std::cout << "\n";
    }
}

//Prints Vector(darray) in a readable format
void printVector(darray A){
    std::vector<double> temp = A.getVector();
    for(int i=0;i<temp.size();i++){
        std::cout << temp[i] << " ";
    }
}


// int main(){
//     Matrix A(3,3);
//     std::vector<double> row1 = {1, 3, 4};
//     A[0] = row1;
//     row1 = {5, 1, 0};
//     A[1] = row1;
//     row1 = {2, 5, 1};
//     A[2] = row1;
//     printMatrix(A);
//     std::cout << "\n";
//     A = reduce(A);
//     printMatrix(A);
    

// }

// PYBIND11 Module definition
PYBIND11_MODULE(operators, m) {
    py::class_<darray>(m, "darray")
        .def(py::init<std::vector<double>>())
        .def("__getitem__", py::overload_cast<int>(&darray::operator[]))  // For getting values
        .def("__setitem__", [](darray& self, int index, double value) { self[index] = value; }) // For setting values
        .def("getSize", &darray::getSize)
        .def("getVector", &darray::getSize) // If you want to expose the underlying vector
        .def("print", &printVector);

    py::class_<Matrix>(m, "Matrix")
        .def(py::init<int, int>())
        .def("__getitem__", py::overload_cast<int>(&Matrix::operator[]))  // For getting rows
        .def("__setitem__", [](Matrix& self, int row, darray value) {
            self[row] = value; 
        }) // For setting rows with a darray
        .def("print", &printMatrix)
        .def("rows", &Matrix::rows)  // Get the number of rows
        .def("cols", &Matrix::cols);  // Get the number of columns
    // Register operations for matrices and darray
    m.def("add", py::overload_cast<darray, darray>(&add));
    m.def("subtract", py::overload_cast<darray, darray>(&subtract));
    m.def("add", py::overload_cast<Matrix, Matrix>(&add));
    m.def("subtract", py::overload_cast<Matrix, Matrix>(&subtract));
    m.def("multiply", py::overload_cast<double, Matrix>(&multiply));
    m.def("multiply", py::overload_cast<Matrix, Matrix>(&multiply));
    m.def("multiply", py::overload_cast<Matrix, darray>(&multiply));
    m.def("reduce", &reduce);
}