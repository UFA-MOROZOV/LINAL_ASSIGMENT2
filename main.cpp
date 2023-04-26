#include "vector"
#include "iostream"
#include "iomanip"
#include "math.h"

using namespace std;


template<typename T>
class Matrix {
public:
    int height = 0;
    int width = 0;
    vector<vector<T>> matrix;

    Matrix() {
        matrix = *new vector<vector<T>>;
    }

    Matrix(int n1, int m1) {
        height = n1;
        width = m1;
        matrix.resize(height);
        for (int i = 0; i < height; i++) {
            matrix[i].resize(width);
        }
    }


    Matrix &operator=(Matrix matrix2) {
        this->matrix.resize(matrix2.height);
        this->height = matrix2.height;
        this->width = matrix2.width;
        for (int i = 0; i < matrix2.height; i++) {
            this->matrix[i].resize(matrix2.width);
            for (int j = 0; j < matrix2.width; j++) {
                this->matrix[i][j] = matrix2.matrix[i][j];
            }
        }
        return *this;
    }


    friend istream &operator>>(istream &input, Matrix &matrix1) {
        int n, m;
        if (matrix1.height == 0 && matrix1.width == 0) {
            input >> n >> m;
            matrix1 = Matrix(n, m);
        }
        for (int i = 0; i < matrix1.height; i++) {
            for (int j = 0; j < matrix1.width; j++) {
                input >> matrix1.matrix[i][j];
            }
        }
        return input;
    }

    friend ostream &operator<<(ostream &output, Matrix matrix1) {
        for (int i = 0; i < matrix1.height; i++) {
            for (int j = 0; j < matrix1.width; j++) {
                if (abs(matrix1.matrix[i][j]) < 1e-10) {
                    output << (T) 0 << " ";
                } else {
                    output << matrix1.matrix[i][j] << " ";
                }
            }
            output << "\n";
        }
        return output;
    }

    friend Matrix operator+(Matrix matrix1, Matrix matrix2) {
        if ((matrix1.height == matrix2.height && matrix1.width == matrix2.width)) {
            Matrix sum = Matrix(matrix1.height, matrix1.width);
            for (int i = 0; i < matrix1.height; i++) {
                for (int j = 0; j < matrix1.width; j++) {
                    sum.matrix[i][j] = matrix1.matrix[i][j] + matrix2.matrix[i][j];
                }
            }
            return sum;
        }
        cout << "Error: the dimensional problem occurred\n";
        return Matrix();
    }

    friend Matrix operator-(Matrix matrix1, Matrix matrix2) {
        if ((matrix1.height == matrix2.height && matrix1.width == matrix2.width)) {
            Matrix dif = Matrix(matrix1.height, matrix1.width);
            for (int i = 0; i < matrix1.height; i++) {
                for (int j = 0; j < matrix1.width; j++) {
                    dif.matrix[i][j] = matrix1.matrix[i][j] - matrix2.matrix[i][j];
                }
            }
            return dif;
        }
        cout << "Error: the dimensional problem occurred\n";
        return Matrix();
    }


    friend Matrix operator*(Matrix matrix1, Matrix matrix2) {
        if (matrix1.width == matrix2.height) {
            Matrix product = Matrix(matrix1.height, matrix2.width);
            T temp;
            for (int i = 0; i < matrix1.height; i++) {
                for (int j = 0; j < matrix2.width; j++) {
                    temp = 0;
                    for (int k = 0; k < matrix1.width; k++) {
                        temp += matrix1.matrix[i][k] * matrix2.matrix[k][j];
                    }
                    product.matrix[i][j] = temp;
                }
            }
            return product;
        }
        cout << "Error: the dimensional problem occurred\n";
        return Matrix();
    }


    Matrix transposed() {
        Matrix answer = Matrix(width, height);
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                answer.matrix[i][j] = matrix[j][i];
            }
        }
        return answer;
    }


    int indexPivot(int column) {
        int maxIndex = column;
        T max1 = matrix[column][column];
        for (int i = column; i < height; i++) {
            if (abs(matrix[i][column]) > max1) {
                maxIndex = i;
                max1 = abs(matrix[i][column]);
            }
        }
        return maxIndex;
    }

};


template<typename T>
class EliminationMatrix;

template<typename T>
class PermutationMatrix;

template<typename T>
class AugmentedMatrix;

template<typename T>
class NormalizationMatrix;


template<typename T>
class SquareMatrix : public Matrix<T> {
public:
    SquareMatrix(Matrix<T> matrix1) {
        this->matrix = matrix1.matrix;
        this->height = matrix1.height;
        this->width = matrix1.width;
    }

    SquareMatrix() : Matrix<T>() {}

    SquareMatrix(int n1) : Matrix<T>(n1, n1) {}

    friend istream &operator>>(istream &input, SquareMatrix &matrix1) {
        int n;
        input >> n;
        matrix1 = SquareMatrix(n);
        for (int i = 0; i < matrix1.height; i++) {
            for (int j = 0; j < matrix1.width; j++) {
                input >> matrix1.matrix[i][j];
            }
        }
        return input;
    }

    friend ostream &operator<<(ostream &output, SquareMatrix matrix1) {
        for (int i = 0; i < matrix1.height; i++) {
            for (int j = 0; j < matrix1.width; j++) {
                if (abs(matrix1.matrix[i][j]) < 1e-10) {
                    output << (T) 0 << " ";
                } else {
                    output << matrix1.matrix[i][j] << " ";
                }
            }
            output << "\n";
        }
        return output;
    }

    friend SquareMatrix operator+(SquareMatrix matrix1, SquareMatrix matrix2) {
        Matrix<T> *m1 = &matrix1;
        Matrix<T> *m2 = &matrix2;
        Matrix sum = *m1 + *m2;
        auto *matrix = static_cast<SquareMatrix *>(&sum);
        return *matrix;
    }

    friend SquareMatrix operator-(SquareMatrix matrix1, SquareMatrix matrix2) {
        Matrix<T> *m1 = &matrix1;
        Matrix<T> *m2 = &matrix2;
        Matrix difference = *m1 - *m2;
        auto *matrix = static_cast<SquareMatrix *>(&difference);
        return *matrix;
    }

    friend SquareMatrix operator*(SquareMatrix matrix1, SquareMatrix matrix2) {
        Matrix<T> *m1 = &matrix1;
        Matrix<T> *m2 = &matrix2;
        Matrix product = (*m1) * (*m2);
        auto *matrix = static_cast<SquareMatrix *>(&product);
        return *matrix;
    }

    SquareMatrix transposed() {
        Matrix<T> *m1 = this;
        Matrix transposed = m1->transposed();
        auto *matrix = static_cast<SquareMatrix *>(&transposed);
        return *matrix;
    }
};


template<typename T>
class IdentityMatrix : public SquareMatrix<T> {
public:
    IdentityMatrix() : SquareMatrix<T>() {}

    explicit IdentityMatrix(int n1) : SquareMatrix<T>(n1) {
        for (int i = 0; i < this->height; i++) {
            this->matrix[i][i] = 1;
        }
    }
};


template<typename T>
class EliminationMatrix : public SquareMatrix<T> {
public:
    EliminationMatrix() : SquareMatrix<T>() {}

    EliminationMatrix(int n1, int n2, SquareMatrix<T> matrix1) : SquareMatrix<T>(matrix1.height) {
        for (int i = 0; i < this->height; i++) {
            this->matrix[i][i] = 1;
        }
        this->matrix[n1][n2] = -matrix1.matrix[n1][n2] / matrix1.matrix[n2][n2];
    }
};


template<typename T>
class PermutationMatrix : public SquareMatrix<T> {
public:
    PermutationMatrix() : SquareMatrix<T>() {}

    PermutationMatrix(int n1, int n2, SquareMatrix<T> matrix1) : SquareMatrix<T>(matrix1.height) {
        for (int i = 0; i < this->height; i++) {
            if (i != n1 && i != n2) {
                this->matrix[i][i] = 1;
            } else if (i == n1) {
                this->matrix[i][n2] = 1;
            } else if (i == n2) {
                this->matrix[i][n1] = 1;
            }
        }
    }
};


template<typename T>
class ColumnVector : public Matrix<T> {
public:
    ColumnVector() = default;

    explicit ColumnVector(int n) {
        this->height = n;
        this->width = 1;
        this->matrix.resize(this->height);
        for (int i = 0; i < this->height; i++) {
            this->matrix[i].resize(this->width);
        }
    };

    explicit ColumnVector(Matrix<T> matrix1) {
        this->height = matrix1.height;
        this->width = matrix1.width;
        this->matrix = matrix1.matrix;
    }

    friend istream &operator>>(istream &input, ColumnVector &matrix1) {
        int n;
        input >> n;
        matrix1 = ColumnVector(n);
        for (int i = 0; i < matrix1.height; i++) {
            input >> matrix1.matrix[i][0];
        }
        return input;
    }

    T norm() {
        T norm = 0;
        for (int i = 0; i < this->height; i++) {
            norm += this->matrix[i][0] * this->matrix[i][0];
        }
        return norm;
    }

    friend ColumnVector operator+(ColumnVector matrix1, ColumnVector matrix2) {
        Matrix<T> *m1 = &matrix1;
        Matrix<T> *m2 = &matrix2;
        Matrix sum = *m1 + *m2;
        auto *matrix = static_cast<ColumnVector *>(&sum);
        return *matrix;
    }

    friend ColumnVector operator-(ColumnVector matrix1, ColumnVector matrix2) {
        Matrix<T> *m1 = &matrix1;
        Matrix<T> *m2 = &matrix2;
        Matrix difference = *m1 - *m2;
        auto *matrix = static_cast<ColumnVector *>(&difference);
        return *matrix;
    }

    friend ColumnVector operator*(Matrix<T> matrix1, ColumnVector matrix2) {
        Matrix<T> *m1 = &matrix1;
        Matrix<T> *m2 = &matrix2;
        Matrix difference = *m1 * *m2;
        auto *matrix = static_cast<ColumnVector *>(&difference);
        return *matrix;
    }

    using Matrix<T>::operator=;
};


template<typename T>
class AugmentedMatrix {
public:
    SquareMatrix<T> A;
    Matrix<T> B;

    AugmentedMatrix(SquareMatrix<T> augmentedMatrix1, SquareMatrix<T> augmentedMatrix2) {
        A = augmentedMatrix1;
        B = augmentedMatrix2;
    }

    friend ostream &operator<<(ostream &output, AugmentedMatrix augmented) {
        for (int i = 0; i < augmented.A.height; i++) {
            for (int j = 0; j < augmented.A.width; j++) {
                output << augmented.A.matrix[i][j] << " ";
            }
            for (int j = 0; j < augmented.B.width; j++) {
                output << augmented.B.matrix[i][j] << " ";
            }
            output << "\n";
        }
        return output;
    }

    friend AugmentedMatrix operator*(Matrix<T> matrix2, AugmentedMatrix matrix1) {
        matrix1.A = SquareMatrix(matrix2) * matrix1.A;
        matrix1.B = matrix2 * matrix1.B;
        return matrix1;
    }

    SquareMatrix<T> inverseMatrix() {
        int step = 1;
        for (int i = 0; i < this->A.width; i++) {
            int pivot = this->A.indexPivot(i);
            if (this->A.indexPivot(i) != i) {
                *this = PermutationMatrix<T>(pivot, i, this->A) * (*this);
                step++;
            }
            for (int j = i + 1; j < this->A.height; j++) {
                if (this->A.matrix[j][i] != 0) {
                    *this = EliminationMatrix<T>(j, i, this->A) * (*this);
                    step++;
                }
            }
        }
        for (int i = this->A.width - 1; i > 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                if (this->A.matrix[j][i] != 0) {
                    *this = EliminationMatrix<T>(j, i, this->A) * *this;
                    step++;
                }
            }
        }
        *this = NormalizationMatrix<T>(this->A) * (*this);
        return B;
    }
};


template<typename T>
class NormalizationMatrix : public SquareMatrix<T> {
public:
    NormalizationMatrix() : SquareMatrix<T>() {}

    explicit NormalizationMatrix(Matrix<T> matrix1) : SquareMatrix<T>(matrix1.height) {
        for (int i = 0; i < this->height; i++) {
            this->matrix[i][i] = (1 / matrix1.matrix[i][i]);
        }
    }

    using Matrix<T>::operator=;
};


template<typename T>
class LeastSquareApproximator {
public:
    void approximate(Matrix<T> coordinates, int step) {
        ColumnVector<double> b = ColumnVector<double>(coordinates.height);
        Matrix<double> A = Matrix<double>(coordinates.height, step + 1);
        for (int i = 0; i < coordinates.height; i++) {
            b.matrix[i][0] = coordinates.matrix[i][1];
            A.matrix[i][0] = 1;
            for (int j = 1; j < step + 1; j++) {
                A.matrix[i][j] = pow(coordinates.matrix[i][0], j);
            }
        }
        cout << "A: " << endl;
        cout << A;
        cout << "A_T*A: " << endl;
        cout << A.transposed() * A;
        cout << "(A_T*A)^-1:" << endl;
        AugmentedMatrix aug = AugmentedMatrix<double>(A.transposed() * A,
                                                      IdentityMatrix<double>(A.width));
        cout << aug.inverseMatrix();
        cout << "A_T*b:" << endl;
        cout << A.transposed() * b;
        cout << "x~:" << endl;
        cout << (Matrix<double>) aug.inverseMatrix() * A.transposed() * b;
    }
};

int main() {
    int n, step;
    cout << fixed << setprecision(4);
    cin >> n;
    Matrix<double> coordinates = Matrix<double>(n, 2);
    cin >> coordinates;
    cin >> step;
    LeastSquareApproximator<double> tool;
    tool.approximate(coordinates, step);
}