#include "s21_matrix_oop.h"

// Private Methods
void S21Matrix::AllocateMemory() {
  matrix_ = new double*[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_]();
  }
}

void S21Matrix::FreeMemory() {
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}

void S21Matrix::resize(int new_rows, int new_cols) {
  // Create a new matrix with the new size
  double** new_matrix = new double*[new_rows];
  for (int i = 0; i < new_rows; i++) {
    new_matrix[i] = new double[new_cols]();
  }

  // Copy the existing data to the new matrix
  int min_rows = std::min(rows_, new_rows);
  int min_cols = std::min(cols_, new_cols);
  for (int i = 0; i < min_rows; i++) {
    for (int j = 0; j < min_cols; j++) {
      new_matrix[i][j] = matrix_[i][j];
    }
  }

  // Delete the old matrix
  for (int i = 0; i < rows_; i++) {
    delete[] matrix_[i];
  }
  delete[] matrix_;

  // Update the matrix pointer and dimensions
  matrix_ = new_matrix;
  rows_ = new_rows;
  cols_ = new_cols;
}

double S21Matrix::GetDeterminant(const S21Matrix& matrix) const {
  double det = 0.0;
  int degree = 1;  // (-1)^(1+j)

  if (matrix.rows_ == 1) {
    det = matrix.matrix_[0][0];
  } else if (matrix.rows_ == 2) {
    det = matrix.matrix_[0][0] * matrix.matrix_[1][1] -
          matrix.matrix_[0][1] * matrix.matrix_[1][0];
  } else {
    for (int j = 0; j < matrix.cols_; j++) {
      S21Matrix tmp = matrix.GetMatrixWithoutRowAndCol(0, j);
      det += (degree * matrix.matrix_[0][j] * GetDeterminant(tmp));
      degree = -degree;
    }
  }
  return det;
}

S21Matrix S21Matrix::GetMatrixWithoutRowAndCol(int row, int col) const {
  S21Matrix result(rows_ - 1, cols_ - 1);
  int offset_row = 0;

  for (int i = 0; i < rows_ - 1; i++) {
    if (i == row) {
      offset_row = 1;
    }
    int offset_col = 0;
    for (int j = 0; j < cols_ - 1; j++) {
      if (j == col) {
        offset_col = 1;
      }
      result.matrix_[i][j] = matrix_[i + offset_row][j + offset_col];
    }
  }
  return result;
}

// Constructors and Destructor
S21Matrix::S21Matrix() : rows_(3), cols_(3) { AllocateMemory(); }

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  AllocateMemory();
}

S21Matrix::S21Matrix(const S21Matrix& other)
    : rows_(other.rows_), cols_(other.cols_) {
  AllocateMemory();
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_), matrix_(other.matrix_) {
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() { FreeMemory(); }

// Accessors and Mutators
int S21Matrix::getRows() const { return rows_; }

int S21Matrix::getCols() const { return cols_; }

void S21Matrix::setRows(int rows) {
  if (rows != rows_) {
    // Handle resizing logic
    resize(rows, cols_);
  }
}

void S21Matrix::setCols(int cols) {
  if (cols != cols_) {
    // Handle resizing logic
    resize(rows_, cols);
  }
}

// Methods
bool S21Matrix::EqMatrix(const S21Matrix& other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    return false;
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      if (std::fabs(matrix_[i][j] - other.matrix_[i][j]) >= 1e-7) {
        return false;
      }
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument(
        "Matrices must have the same dimensions for addition.");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument(
        "Matrices must have the same dimensions for subtraction.");
  }
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument(
        "Number of columns of the first matrix must be equal to the number of "
        "rows of the second matrix.");
  }

  S21Matrix result(rows_, other.cols_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      for (int k = 0; k < cols_; ++k) {
        result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = std::move(result);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::invalid_argument(
        "Matrix must be square to calculate determinant.");
  }
  return GetDeterminant(*this);
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_) {
    throw std::invalid_argument(
        "Matrix must be square to calculate complements.");
  }

  S21Matrix result(rows_, cols_);

  if (rows_ == 1) {
    result.matrix_[0][0] = 1.0;
    return result;
  }

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor = GetMatrixWithoutRowAndCol(i, j);
      double det = GetDeterminant(minor);
      result.matrix_[i][j] = pow(-1, i + j) * det;
    }
  }

  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (rows_ != cols_) {
    throw std::invalid_argument("Matrix must be square to calculate inverse.");
  }

  double det = Determinant();

  if (fabs(det) < 1e-7) {
    throw std::runtime_error(
        "Matrix determinant is zero, inverse cannot be calculated.");
  }

  S21Matrix complements = CalcComplements();
  S21Matrix transposed = complements.Transpose();
  transposed.MulNumber(1.0 / det);
  return transposed;
}

// Operator Overloading
S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result = *this;
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result = *this;
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result = *this;
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(double num) {
  S21Matrix result = *this;
  result.MulNumber(num);
  return result;
}

bool S21Matrix::operator==(const S21Matrix& other) { return EqMatrix(other); }

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (this != &other) {
    FreeMemory();
    rows_ = other.rows_;
    cols_ = other.cols_;
    AllocateMemory();
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(double num) {
  MulNumber(num);
  return *this;
}

double& S21Matrix::operator()(int i, int j) {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::out_of_range("Index is out of range");
  }
  return matrix_[i][j];
}

const double& S21Matrix::operator()(int i, int j) const {
  if (i >= rows_ || j >= cols_ || i < 0 || j < 0) {
    throw std::out_of_range("Index is out of range");
  }
  return matrix_[i][j];
}