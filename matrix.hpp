#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <iostream>
#include <cstring>
#include <chrono>
#include <fstream>
#include <sstream>
#include <array>
#include </opt/homebrew/Cellar/openblas/0.3.23/include/cblas.h>
#include </opt/homebrew/Cellar/libomp/16.0.3/include/omp.h>


template<class T>
class Matrix {
public:
    size_t row;
    size_t col;
    size_t channel;
    size_t spans;
    size_t total;
    T *elements;

    /*>>> static counter >>>*/
    int *ref_cnt;
    static int mat_num;

public:
    /*>>> getters >>>*/
    size_t getRow() {
        return this->row;
    }

    size_t getCol() {
        return this->col;
    }

    size_t getChannel() {
        return this->channel;
    }

    size_t getSpans() {
        return this->spans;
    }

    size_t getTotal() {
        return this->total;
    }

    size_t getSize() {
        return this->row * this->col;
    }

    T getElementValue(int x, int y, int z) const {
        if (x < 0 || y < 0 || z < 0 || x >= this->row || y >= this->col || z >= this->channel) {
            std::cerr << "error: invalid param" << std::endl;
            exit(EXIT_FAILURE);
        }
        return this->elements[x * spans + y * channel + z];
    }

    T &getElementPointer(int x, int y, int z) const {
        if (x < 0 || y < 0 || z < 0 || x >= this->row || y >= this->col || z >= this->channel) {
            std::cerr << "error: invalid param" << std::endl;
            exit(EXIT_FAILURE);
        }
        return this->elements[x * spans + y * channel + z];
    }


    /*>>> setters >>>*/
    void setElements(T *contents) {
//        std::cout << "info: size of construct array " << len << std::endl;
//        std::cout << this->total << std::endl;
//        if (this->total != len) {
//            std::cerr << "error: invalid elements length" << std::endl;
//        }
        int len = this->row * this->col * this->channel;
        for (int i = 0; i <= len; i++) {
            this->elements[i] = contents[i];
        }
    }


    /*>>> static counter, constructor and destructor >>>*/
    int getMatrixNum() {
        return mat_num;
    }

//    Matrix<T> createMatrix(int row = 128, int col = 128, int channel = 3) {
//        std::cout << "info: manuel constructor is invoked" << std::endl;
//        if (row <= 0 || col <= 0 || channel <= 0) {
//            std::cerr << "error: invalid param" << std::endl;
//            exit(EXIT_FAILURE);
//        }
//        Matrix<T> matrix;
//        matrix.row = row;
//        matrix.col = col;
//        matrix.channel = channel;
//        matrix.spans = col * channel;
//        matrix.total = row * col * channel;
//        matrix.elements = new T[row * col * channel];
//        memset(this->elements, 0, sizeof(T) * row * col * channel);
//        mat_num += 1;
//        std::cout << "info: exist " << mat_num << " matrices" << std::endl;
//        return matrix;
//    }

    explicit Matrix(int row = 128, int col = 128, int channel = 3) {
        std::cout << "info: numerical constructor is invoked" << std::endl;
        if (row <= 0 || col <= 0 || channel <= 0) {
            std::cerr << "error: invalid param" << std::endl;
            exit(EXIT_FAILURE);
        }
        this->ref_cnt = new int[1];
        this->ref_cnt[0] = 1;
        this->row = row;
        this->col = col;
        this->channel = channel;
        this->spans = this->col * this->channel;
        this->total = this->row * this->col * this->channel;
        this->elements = new T[row * col * channel];
        memset(this->elements, 0, sizeof(T) * row * col * channel);
        mat_num += 1;
        std::cout << "info: exist " << mat_num << " matrices" << std::endl;
    }

    explicit Matrix(std::ifstream &fin) {
        std::cout << "info: ifstream constructor is invoked" << std::endl;

        std::cout << "info: exist " << mat_num << " matrices" << std::endl;
    }

    Matrix(const Matrix<T> &matrix) {
        std::cout << "info: copy constructor is invoked" << std::endl;
        this->row = matrix.row;
        this->col = matrix.col;
        this->channel = matrix.channel;
        this->ref_cnt = matrix.ref_cnt;
        this->spans = matrix.spans;
        this->elements = matrix.elements;
        mat_num += 1;
        std::cout << "info: exist " << mat_num << " matrices" << std::endl;
    }

    ~Matrix() {
        std::cout << "info: default destructor is invoked" << std::endl;
        mat_num -= 1;
        if (this->ref_cnt[0] == 1) {
            delete[] this->ref_cnt;
            delete[] this->elements;
        } else {
            this->ref_cnt[0] -= 1;
        }
        std::cout << "info: exist " << mat_num << " matrices" << std::endl;
    }


    /*>>> overload operators >>>*/
    Matrix<T> Transpose() {
        Matrix<T> matrix;
        matrix.row = this->row;
        matrix.col = this->col;
        matrix.channel = this->channel;
        matrix.spans = this->spans;
        matrix.total = this->total;
        matrix.elements = new T[this->row * this->col * this->channel];
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    *matrix.getElementPointer(i, j, k) =
                            this->getElementPointer(j, i, k);
                }
            }
        }
        return matrix;
    }

    Matrix<T> operator+(T num) const {
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    this->getElementPointer(i, j, k) += num;
                }
            }
        }
        return *this;
    }

    Matrix<T> operator+(const Matrix<T> &matrix) {
        if (this->row != matrix.row || this->col != matrix.col || this->channel != matrix.channel) {
            std::cerr << "error: invalid param" << std::endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    *this->getElementPointer(i, j, k) +=
                            matrix.getElementValue(i, j, k);
                }
            }
        }
        return *this;
    }

    Matrix<T> operator++() {

    }

    Matrix<T> operator-(T num) const {
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    this->getElementPointer(i, j, k) -= num;
                }
            }
        }
        return *this;
    }

    Matrix<T> operator-(const Matrix<T> &matrix) {
        if (this->row != matrix.row || this->col != matrix.col || this->channel != matrix.channel) {
            std::cerr << "error: invalid param" << std::endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    this->getElementPointer(i, j, k) -=
                            matrix.getElementValue(i, j, k);
                }
            }
        }
        return *this;
    }

    Matrix<T> operator*(T num) const {
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    this->getElementPointer(i, j, k) *= num;
                }
            }
        }
        return *this;
    }

    Matrix<T> operator*(const Matrix<T> &matrix) {
        if (this->row != matrix.row || this->col != matrix.col || this->channel != matrix.channel) {
            std::cerr << "error: invalid param" << std::endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    this->getElementPointer(i, j, k) *=
                            matrix.getElementValue(i, j, k);
                }
            }
        }
        return *this;
    }

    Matrix<T> operator/(T num) const {
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    this->getElementPointer(i, j, k) /= num;
                }
            }
        }
        return *this;
    }

    Matrix<T> operator/(const Matrix<T> &matrix) {
        if (this->row != matrix.row || this->col != matrix.col || this->channel != matrix.channel) {
            std::cerr << "error: invalid param" << std::endl;
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < this->row; i++) {
            for (int j = 0; j < this->col; j++) {
                for (int k = 0; k < this->channel; k++) {
                    this->getElementPointer(i, j, k) /=
                            matrix.getElementValue(i, j, k);
                }
            }
        }
        return *this;
    }
};

template<typename T>
int Matrix<T>::mat_num = 0;

template<typename T>
std::ostream &operator<<(std::ostream &ostream, const Matrix<T> &matrix) {
    ostream << "print: ";
    ostream << "[";
    ostream << "\n";
    for (int k = 0; k < matrix.channel; k++) {
        ostream << "[";
        for (int i = 0; i < matrix.row; i++) {
            ostream << "[";
            for (int j = 0; j < matrix.col; j++) {
                if (j == matrix.col - 1) {
                    ostream << matrix.getElementValue(i, j, k);
                } else {
                    ostream << matrix.getElementValue(i, j, k) << " ";
                }
            }
            ostream << "]";
        }
        ostream << "]";
        ostream << "\n";
    }
    ostream << "]";
    return ostream;
}

template<typename T>
void multiply(const Matrix<T> &matrix1, const Matrix<T> &matrix2, Matrix<T> &result) {
#pragma omp sections
    if (matrix1.row != matrix2.col || matrix1.channel != matrix2.channel
        || result.row != matrix1.row || result.col != matrix2.col || matrix1.channel != result.channel) {
        std::cerr << "error: invalid matrix size" << std::endl;
        exit(EXIT_FAILURE);
    }

//    float *in1 = new float[matrix1.row * matrix1.col * matrix1.channel];
//    memset(in1, 0, sizeof(float) * (matrix1.row * matrix1.col * matrix1.channel));
//    float *in2 = new float[matrix2.row * matrix2.col * matrix2.channel];
//    memset(in2, 0, sizeof(float) * (matrix2.row * matrix2.col * matrix2.channel));
    T *res = new T[result.row * result.col * result.channel];
    memset(res, 0, sizeof(T) * (result.row * result.col * result.channel));

    int cnt1 = 0;
    int cnt2 = 0;
    int cnt3 = 0;

//    for (int k = 0; k < matrix1.channel; k++) {
//        for (int i = 0; i < matrix1.row; i++) {
//            for (int j = 0; j < matrix1.col; j++) {
//                std::cout << res[cnt1] << " ";
//                std::cout << matrix1.getElementValue(i, j, k) << " ";
//                in1[cnt1++] =
//                        matrix1.getElementValue(i, j, k);
//            }
//        }
//    }
//    std::cout << "\n";
//    for (int k = 0; k < matrix2.channel; k++) {
//        for (int i = 0; i < matrix2.row; i++) {
//            for (int j = 0; j < matrix2.col; j++) {
//                std::cout << res[cnt2] << " ";
//                std::cout << matrix2.getElementValue(i, j, k) << " ";
//                in2[cnt2++] =
//                        matrix2.getElementValue(i, j, k);
//            }
//        }
//    }
//    std::cout << "\n";

    int row = matrix1.row;
    int col = matrix1.col;

//    std::cout << row << std::endl;
//    std::cout << col << std::endl;

    if constexpr (std::is_same<T, int>::value) {
        cblas_sgemm(
                CblasRowMajor, CblasNoTrans, CblasNoTrans,
                matrix1.row, matrix2.col, matrix1.col, 1.0f,
                (float *)matrix1.elements, matrix1.col,
                (float *)matrix2.elements, matrix2.col, 0,
                (float *)res, matrix2.col);
    } else if constexpr (std::is_same<T, float>::value) {
        cblas_sgemm(
                CblasRowMajor, CblasNoTrans, CblasNoTrans,
                matrix1.row, matrix2.col, matrix1.col, 1.0f,
                matrix1.elements, matrix1.col,
                matrix2.elements, matrix2.col, 0,
                res, matrix2.col);
    } else if constexpr (std::is_same<T, double>::value) {
        cblas_dgemm(
                CblasRowMajor, CblasNoTrans, CblasNoTrans,
                matrix1.row, matrix2.col, matrix1.col, 1.0,
                matrix1.elements, matrix1.col,
                matrix2.elements, matrix2.col, 0,
                res, matrix2.col);
    } else {
        std::cerr << "error: unsupported data type" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < result.row; i++) {
        for (int j = 0; j < result.col; j++) {
            for (int k = 0; k < result.channel; k++) {
//                std::cout << cnt3 << ":" << res[cnt3] << " ";
                result.getElementPointer(i, j, k) = res[cnt3++];
            }
        }
    }
//    std::cout << "\n";
//    delete[] in1;
//    delete[] in2;
    delete[] res;
}

#endif