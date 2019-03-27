#ifndef __MatriX__
#define __MatriX__

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <stdint.h>
#include <vector>

namespace mx {

struct Dimensions {
    uint32_t rows;
    uint32_t collumns;
};

bool operator==(const Dimensions &a, const Dimensions &b) {
    if (a.rows == b.rows && a.collumns == b.collumns)
        return true;
    else
        return false;
}

template <typename T>
class Matrix {
  public:
  
  /**
  * To create a matrix, we need to know its size
  * Therefore the default constructor is disabled
  **/
    Matrix() = delete;

  /**
  * Constructs a matrix of given size 
  * and sets its elements to default value of their type
  **/
    Matrix(Dimensions d) : dim(d) {
        elements.resize(dim.rows * dim.collumns, T());
    }
    
  /**
  * Constructs a matrix of given size 
  * and sets its elements to values in initializer list.
  * (Starting from first row)
  * Throws an exception if the number of elements in the list
  * Is not equal to number of elements in the matrix.
  **/
    Matrix(Dimensions d, std::initializer_list<T> l) : dim(d), elements(l) {
        if ((dim.rows * dim.collumns) != l.size())
            throw std::invalid_argument("Matrix: Make sure that row_count * "
                                        "collumn_count = list.lenght()");
    }
    
  /**
  * Constructs a matrix with one row and fills it with 
  * values from list l
  **/
    Matrix(std::initializer_list<T> l) : elements(l) {
        dim.rows = 1;
        dim.collumns = l.size();
    }

    Matrix(const Matrix &m) = default;
    Matrix &operator=(const Matrix &m) = default;
    Matrix(Matrix &&m) = default;
    Matrix &operator=(Matrix &&m) = default;
    ~Matrix() = default;

    bool isSquare() const { return (dim.collumns == dim.rows); }

 /** 
 * Returns reference to an element ofbthe matrix.
 * Indexes start at 1
 **/
    T &at(uint32_t row, uint32_t collumn) {
        if (row > dim.rows || collumn > dim.collumns || row < 1 ||
            collumn < 1) {
            throw std::out_of_range("Matrix.at() exception: Invalid index");
        }
        try {
            return elements.at((row - 1) * dim.collumns + (collumn - 1));
        } catch (std::out_of_range e) {
            std::cerr << "Matrix.at() Out of range Exception: " << e.what()
                      << "\n";
            throw;
        } catch (std::exception e) {
            std::cerr << "Unknown Exception in Matrix.at(): " << e.what()
                      << "\n";
            throw;
        }
    }

    Matrix subMatrix(Dimensions start, Dimensions range) {
        if ((start.collumns + range.collumns) > dim.collumns ||
            (start.rows + range.rows) > dim.rows || start.collumns < 1 ||
            start.rows < 1)
            throw std::invalid_argument("Cannot create submatrix of this size");

        Matrix nMatrix{range};
        for (int i = start.rows; i <= range.rows; i++) {
            for (int j = start.collumns; j <= range.collumns; j++) {
                nMatrix.at(i, j) = this->at(i, j);
            }
        }

        return nMatrix;
    }

    Matrix minor(Dimensions number) {
        if (number.collumns > dim.collumns || number.rows > dim.rows ||
            number.collumns < 1 || number.rows < 1)
            throw std::invalid_argument("Cannot create minor!");

        Dimensions mDim{dim.rows - 1, dim.collumns - 1};
        Matrix minor{mDim};
        int k = 1;

        for (int i = 1; i <= dim.rows; i++) {
            int l = 1;
            if (i == number.rows)
                continue;

            for (int j = 1; j <= dim.collumns; j++) {
                if (j == number.collumns)
                    continue;

                minor.at(k, l) = this->at(i, j);
                l++;
            }
            k++;
        }
        return minor;
    }
    
    /**
    * Calculates the determinant recursively using
    * Laplace expansion
    **/

    T det() {
        if (!isSquare())
            throw std::invalid_argument(
                "Cannot calculate determinant of a non-square matrix");

        if (dim.collumns == 1) {
            return at(1, 1);
        }

        T sum = T();

        for (uint32_t j = 1; j <= dim.collumns; j++) {

            sum += at(1, j) * pow(-1, 1 + j) * minor({1, j}).det();
        }

        return sum;
    }
    
    /** Returns a new matrix that 
    * is a transposed parent matrix
    **/
    
    Matrix transpose() {
        Matrix transposed{dim};

        for (uint32_t i = 1; i <= dim.rows; i++) {
            for (uint32_t j = 1; j <= dim.collumns; j++) {
                transposed.at(j, i) = at(i, j);
            }
        }

        return transposed;
    }
    
    /** Returns a new matrix that 
    * is an inverted parent matrix
    **/
    Matrix invert() {
        if (!isSquare())
            throw std::invalid_argument("Matrix must be square to be inverted");

        T determinant = det();

        if (determinant == 0)
            throw std::invalid_argument(
                "Cannot invert matrix, determinant is equal to zero");

        Matrix inverted{dim};

        for (uint32_t i = 1; i <= dim.rows; i++) {
            for (uint32_t j = 1; j <= dim.collumns; j++) {
                inverted.at(j, i) =
                    minor({i, j}).det() * (1 / determinant) * pow(-1, i + j);
            }
        }

        return inverted;
    }

    const Dimensions &getDimensions() const { return dim; }

    const std::vector<T> &getVector() const { return elements; }

  private:
    Dimensions dim;
    std::vector<T> elements;
};

//Arithmetic operators overloads

template <typename T>
Matrix<T> operator+(Matrix<T> &a, Matrix<T> &b) {
    if (a.getDimensions() == b.getDimensions()) {
        auto d = a.getDimensions();
        Matrix<T> mx{d};
        for (int i = 1; i <= d.rows; i++) {
            for (int j = 1; j <= d.collumns; j++) {
                mx.at(i, j) = a.at(i, j) + b.at(i, j);
            }
        }
        return mx;
    }

    throw std::invalid_argument("Attempted to add incompatible matrices");
}

template <typename T>
Matrix<T> operator-(Matrix<T> &a, Matrix<T> &b) {
    return a + (b * -1.0);
}

template <typename T>
Matrix<T> operator*(Matrix<T> &a, T scalar) {
    auto d = a.getDimensions();
    Matrix<T> newMatrix = a;
    for (int i = 1; i <= d.rows; i++) {
        for (int j = 1; j <= d.collumns; j++) {
            newMatrix.at(i, j) *= scalar;
        }
    }
    return newMatrix;
}

template <typename T>
Matrix<T> operator*(Matrix<T> &a, Matrix<T> &b) {
    if (a.getDimensions().collumns == b.getDimensions().rows) {
        auto d = b.getDimensions();
        Matrix<T> mx{
            Dimensions{a.getDimensions().rows, b.getDimensions().collumns}};
        for (int i = 1; i <= d.rows; i++) {
            for (int j = 1; j <= d.collumns; j++) {
                mx.at(i, j) = 0;
                for (int k = 1; k <= d.rows; k++) {
                    mx.at(i, j) += a.at(i, k) * b.at(k, j);
                }
            }
        }
        return mx;
    }

    throw std::invalid_argument("Attempted to multiply incompatible matrices");
}
} // namespace mx

#endif