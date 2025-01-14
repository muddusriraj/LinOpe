# LinOpeC++

Simple C++ Library (also compiled into a Python extension using pybind11) for Basic Matrix Operations and Vector plotting.

Includes:

(MatrixOperations/operators.h)
"darray" object type used to store vectors which can hold doubles.
-   Wrapper for std::vector<double> to allow for use in Python
-   Can be initialized by passing an std::vector into the constructor OR a darray size to initialize all values to 0
-   Contains:
    -   [] operator for indexing (can be used to retrieve and set values)
    -   getSize() to retrieve size of darray
    -   \* operator which can be used for darray * (darray) multiplication (with an output of type double) or darray * (double) multiplication with output of type darray.
    -   \+ operator which can be used for darray + (darray) addition (with output of type darray)
    -   \- operator which can be used for darray - (darray) subtraction (with output of type darray)
    -   norm() to retrieve the norm of the vector

"Matrix" object type used to store vectors which hold type darray.
-   Stores std::vector<darray>
-   Can be initialized with the address of another Matrix (copies all the values in) or dimensions eg. Matrix{3, 3} for a 3x3 Matrix
-   Contains:
    -   [] operator to retrieve and set darrays inside matrix
        - individual values can be retrieved / set by chaining these operators eg. x[1][2] = 3;
    -   rows() retrieves number of rows
    -   cols() retrieves number of columns
    -   reduce() transforms Matrix to row echelon form (refer to standalone rowReduce method for RREF)
    -   \* operator for Matrix * (Matrix) -> Matrix, Matrix * (darray) -> darray, Matrix * (double) -> Matrix

Standalone Methods:
-   darray subtract(darray A, darray B)
-   darray add(darray A, darray B)
-   Matrix add(Matrix A, Matrix B)
-   Matrix subtract(Matrix A, Matrix B)
-   darray multiply(double scalar, darray vectA)
-   Matrix multiply(double scalar, Matrix A)
-   Matrix multiply(Matrix A, Matrix B)
-   darray multiply(Matrix A, darray vectA)
-   Matrix reduce(Matrix A)
-   Matrix rowReduce(Matrix A)
-   Matrix transpose(Matrix A)
-   Matrix identity(int x)
-   double det(Matrix A)
-   darray normalize(darray A)
-   void printMatrix(Matrix A)
-   void printVector(darray A)

(MatrixOperations/graph.h)
Methods:
-   void plot(const std::vector<darray>& vectors, const std::vector<std::string> colors = {"blue"})
        - Allows user to plot vectors
