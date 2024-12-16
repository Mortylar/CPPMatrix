#ifndef S21_MATRIX_OOP_H_
#define S21_MATRIX_OOP_H_

class S21Matrix {
 public:
  S21Matrix();

  S21Matrix(int rows, int cols);

  S21Matrix(const S21Matrix& other);

  S21Matrix(S21Matrix&& other);

  S21Matrix& operator=(const S21Matrix& other);

  S21Matrix& operator=(S21Matrix&& other);

  ~S21Matrix();
  void FirstLayerInitialize(int rows, int cols);
  void SecondLayerInitialize(int rows, int cols);
  void S21PutinFunction();

  //=======accessors_and_mutators==========//

  int GetRows();
  int GetCols();

  void SetRows(int rows);
  void SetCols(int cols);

  void S21ResizeMatrix(int rows, int cols);
  //================methods===============//

  bool EqMatrix(const S21Matrix& other);

  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  bool CheckSumMatrixSize(const S21Matrix& other);
  void SumOrNotSumMatrix(const S21Matrix& other, bool is_size);

  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  bool CheckMulMatrixSize(const S21Matrix& other);

  S21Matrix Transpose();

  S21Matrix CalcComplements();
  double Determinant();
  double S21BabyDeterminant();
  S21Matrix S21GetMinorMatrix(int row, int col);
  bool CheckDeterminant();

  S21Matrix InverseMatrix();

  //===============overload================//

  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  S21Matrix operator+=(const S21Matrix& other);
  S21Matrix operator-=(const S21Matrix& other);
  S21Matrix operator*=(const S21Matrix& other);

  bool const operator==(const S21Matrix& other);
  double& operator()(int row, int col);

  friend S21Matrix operator*(const double num, const S21Matrix& other);
  friend S21Matrix operator*(const S21Matrix& other, const double num);

  //===============private_fields==========//
 private:
  int rows_;
  int cols_;
  double* matrix_;
};

S21Matrix operator*(const double num, const S21Matrix& other);
S21Matrix operator*(const S21Matrix& other, const double num);

//==============math=============//

#define S21_ANTI_EPS 10e7

double S21DAbs(double first, double second);
double S21DMin(double first, double second);
bool S21EqDouble(double first, double second);

#endif  //  S21_MATRIX_OOP_H_
