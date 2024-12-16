#include "s21_matrix_oop.h"

#include <cstring>
#include <iostream>

S21Matrix::S21Matrix() { this->SecondLayerInitialize(1, 1); }

S21Matrix::S21Matrix(int rows, int cols) {
  if ((rows <= 0) || (cols <= 0)) {
    throw std::out_of_range("Incorrect matrix size");
  }
  this->SecondLayerInitialize(rows, cols);
}

S21Matrix::S21Matrix(const S21Matrix& other) {
  if (&other != this) {
    this->FirstLayerInitialize(other.rows_, other.cols_);
    std::memcpy(this->matrix_, other.matrix_,
                (this->rows_ * this->cols_ * sizeof(double)));
  }
}

S21Matrix::S21Matrix(S21Matrix&& other) {
  if (&other != this) {
    this->FirstLayerInitialize(other.rows_, other.cols_);
    std::swap(this->matrix_, other.matrix_);
  }
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  if (&other != this) {
    this->~S21Matrix();
    this->FirstLayerInitialize(other.rows_, other.cols_);
    std::memcpy(this->matrix_, other.matrix_,
                (this->rows_ * this->cols_ * sizeof(double)));
  }
  return *this;
}

S21Matrix& S21Matrix::operator=(S21Matrix&& other) {
  if (&other != this) {
    this->~S21Matrix();
    this->FirstLayerInitialize(other.rows_, other.cols_);
    std::swap(this->matrix_, other.matrix_);
  }
  return *this;
}

S21Matrix::~S21Matrix() {
  this->rows_ = 0;
  this->cols_ = 0;
  if (this->matrix_) {
    delete[] this->matrix_;
    this->matrix_ = NULL;
  }
}

void S21Matrix::FirstLayerInitialize(int rows, int cols) {
  this->rows_ = rows;
  this->cols_ = cols;
  this->matrix_ = new double[rows * cols]();
  return;
}

void S21Matrix::SecondLayerInitialize(int rows, int cols) {
  this->FirstLayerInitialize(rows, cols);
  this->S21PutinFunction();
  return;
}

void S21Matrix::S21PutinFunction() {
  for (int i = 0; i < (this->rows_ * this->cols_); ++i) {
    this->matrix_[i] = 0.0;
  }
  return;
}

//===========accessors_and_mutators=========//

int S21Matrix::GetRows() { return (this->rows_); }

int S21Matrix::GetCols() { return (this->cols_); }

void S21Matrix::SetRows(int rows) {
  if (rows <= 0) {
    throw std::out_of_range("Incorrect rows number");
  }
  if (rows != this->rows_) {
    this->S21ResizeMatrix(rows, this->cols_);
  }
  return;
}

void S21Matrix::SetCols(int cols) {
  if (cols <= 0) {
    throw std::out_of_range("Incorrect cols number");
  }
  if (cols != this->cols_) {
    this->S21ResizeMatrix(this->rows_, cols);
  }
  return;
}

void S21Matrix::S21ResizeMatrix(int rows, int cols) {
  S21Matrix tmp_matrix(rows, cols);
  int min_rows = ((rows > this->rows_) ? (this->rows_) : (rows));
  int min_cols = ((cols > this->cols_) ? (this->cols_) : (cols));
  for (int i = 0; i < min_rows; ++i) {
    for (int j = 0; j < min_cols; ++j) {
      tmp_matrix.matrix_[i * tmp_matrix.cols_ + j] =
          this->matrix_[i * this->cols_ + j];
    }
  }
  *this = tmp_matrix;
}

//===================methods===============//

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  bool is_eq = true;
  if ((this->rows_ == other.rows_) && (this->cols_ == other.cols_)) {
    for (int i = 0; (is_eq && (i < (this->rows_ * this->cols_))); ++i) {
      is_eq = S21EqDouble(this->matrix_[i], other.matrix_[i]);
    }
  } else {
    is_eq = false;
  }
  return is_eq;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (!(this->CheckSumMatrixSize(other))) {
    throw std::logic_error("Unequal matrix size");
  }
  bool is_sum = true;
  this->SumOrNotSumMatrix(other, is_sum);
  return;
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (!(this->CheckSumMatrixSize(other))) {
    throw std::logic_error("Unequal matrix size");
  }
  bool is_sum = false;
  this->SumOrNotSumMatrix(other, is_sum);
  return;
}

bool S21Matrix::CheckSumMatrixSize(const S21Matrix& other) {
  return ((this->rows_ == other.rows_) && (this->cols_ == other.cols_));
}

void S21Matrix::SumOrNotSumMatrix(const S21Matrix& other, bool is_sum) {
  for (int i = 0; i < (this->rows_ * this->cols_); ++i) {
    if (is_sum) {
      this->matrix_[i] += other.matrix_[i];
    } else {
      this->matrix_[i] -= other.matrix_[i];
    }
  }
  return;
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < (this->rows_ * this->cols_); ++i) {
    this->matrix_[i] *= num;
  }
  return;
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  if (!CheckMulMatrixSize(other)) {
    throw std::logic_error("Incorrect matrix size");
  }
  S21Matrix result{this->rows_, other.cols_};
  for (int i = 0; i < result.rows_; ++i) {
    for (int j = 0; j < result.cols_; ++j) {
      for (int k = 0; k < this->cols_; ++k) {
        result.matrix_[i * result.cols_ + j] +=
            (this->matrix_[i * this->cols_ + k] *
             other.matrix_[k * other.cols_ + j]);
      }
    }
  }
  *this = result;
  return;
}

bool S21Matrix::CheckMulMatrixSize(const S21Matrix& other) {
  return (this->cols_ == other.rows_);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result{this->cols_, this->rows_};
  for (int i = 0; i < result.rows_; ++i) {
    for (int j = 0; j < result.cols_; ++j) {
      result.matrix_[i * result.cols_ + j] = this->matrix_[j * this->cols_ + i];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix result(1, 1);
  if (!this->CheckDeterminant()) {
    throw std::logic_error("Incorrect matrix size");
  }
  if (this->rows_ == 1) {
    result.matrix_[0] = 1.0;
  } else {
    result = *this;
    for (int i = 0; i < this->rows_; ++i) {
      for (int j = 0; j < this->cols_; ++j) {
        S21Matrix minor_matrix = this->S21GetMinorMatrix(i, j);
        double minor_number = minor_matrix.Determinant();
        int doubting_unit = (((i + j) % 2) ? (-1) : 1);
        result.matrix_[i * result.cols_ + j] = (doubting_unit * minor_number);
      }
    }
  }
  return result;
}

double S21Matrix::Determinant() {
  double s21_determinant = 0.0;
  if (!this->CheckDeterminant()) {
    throw std::logic_error("Incorrect matrix size");
  }
  int baby_size = 3;
  if (this->rows_ <= baby_size) {
    s21_determinant = this->S21BabyDeterminant();
  } else {
    int doubting_unit = 1;
    for (int i = 0; i < this->rows_; ++i) {
      S21Matrix minor_matrix = this->S21GetMinorMatrix(i, 0);
      double local_det = minor_matrix.Determinant();
      s21_determinant +=
          (doubting_unit * local_det * this->matrix_[i * this->cols_]);
      doubting_unit *= (-1);
    }
  }
  return s21_determinant;
}

double S21Matrix::S21BabyDeterminant() {
  double det = 0.0;
  int matrix_size = this->rows_;
  if (1 == matrix_size) {
    det = this->matrix_[0];
  } else if (2 == matrix_size) {
    det = (this->matrix_[0] * this->matrix_[3]);
    det -= (this->matrix_[1] * this->matrix_[2]);
  } else if (3 == matrix_size) {
    det += (this->matrix_[0] * this->matrix_[4] * this->matrix_[8]);
    det += (this->matrix_[1] * this->matrix_[5] * this->matrix_[6]);
    det += (this->matrix_[2] * this->matrix_[3] * this->matrix_[7]);
    det -= (this->matrix_[2] * this->matrix_[4] * this->matrix_[6]);
    det -= (this->matrix_[0] * this->matrix_[5] * this->matrix_[7]);
    det -= (this->matrix_[1] * this->matrix_[3] * this->matrix_[8]);
  }
  return det;
}

S21Matrix S21Matrix::S21GetMinorMatrix(int row, int col) {
  S21Matrix minor_matrix((this->rows_ - 1), (this->cols_ - 1));
  for (int i = 0; i < minor_matrix.rows_; ++i) {
    for (int j = 0; j < minor_matrix.cols_; ++j) {
      if ((i < row) && (j < col)) {
        minor_matrix.matrix_[i * minor_matrix.cols_ + j] =
            this->matrix_[i * this->cols_ + j];
      } else if ((i < row) && (j >= col)) {
        minor_matrix.matrix_[i * minor_matrix.cols_ + j] =
            this->matrix_[i * this->cols_ + j + 1];
      } else if ((i >= row) && (j < col)) {
        minor_matrix.matrix_[i * minor_matrix.cols_ + j] =
            this->matrix_[(i + 1) * this->cols_ + j];
      } else {
        minor_matrix.matrix_[i * minor_matrix.cols_ + j] =
            this->matrix_[(i + 1) * this->cols_ + j + 1];
      }
    }
  }
  return minor_matrix;
}

bool S21Matrix::CheckDeterminant() { return (this->cols_ == this->rows_); }

S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix result;
  double this_det = this->Determinant();
  if (S21EqDouble(this_det, 0.0)) {
    throw std::runtime_error("Division by zero");
  }
  result = this->CalcComplements();
  result = result.Transpose();
  result.MulNumber((1.0 / this_det));
  return result;
}

//===================overload========================//

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix result = (*this);
  result.SumMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix result = (*this);
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix result = (*this);
  result.MulMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator+=(const S21Matrix& other) {
  this->SumMatrix(other);
  return (*this);
}

S21Matrix S21Matrix::operator-=(const S21Matrix& other) {
  this->SubMatrix(other);
  return (*this);
}

S21Matrix S21Matrix::operator*=(const S21Matrix& other) {
  this->MulMatrix(other);
  return (*this);
}

bool const S21Matrix::operator==(const S21Matrix& other) {
  return (this->EqMatrix(other));
}

double& S21Matrix::operator()(int row, int col) {
  if ((row >= this->rows_) || (row < 0) || (col >= this->cols_) || (col < 0)) {
    throw std::out_of_range("Incorrect argument");
  }
  return *(this->matrix_ + (row * this->cols_ + col));
}

S21Matrix operator*(const double num, const S21Matrix& matrix) {
  S21Matrix result = matrix;
  result.MulNumber(num);
  return result;
}

S21Matrix operator*(const S21Matrix& matrix, const double num) {
  S21Matrix result = matrix;
  result.MulNumber(num);
  return result;
}

//=======================math======================//

double S21DAbs(double d_number) {
  return ((d_number >= 0.0) ? (d_number) : (-d_number));
}

double S21DMin(double first, double second) {
  return ((first <= second) ? (first) : (second));
}

bool S21EqDouble(double first, double second) {
  return ((S21DAbs(first - second) * S21_ANTI_EPS)) <=
         S21DMin(S21DAbs(first), S21DAbs(second));
}
