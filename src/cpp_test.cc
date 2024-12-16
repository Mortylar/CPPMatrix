#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(IsEqDoubleTest, EqUnit) {
  double s21_eps = 1e-7;
  double start_x = 1e-4;
  double shifting_koef = 0.09;
  for (int i = 0; i < 100000000; ++i) {
    double delta_x = start_x * 1e-2;
    double local_eps = start_x * s21_eps;
    for (double x = start_x; x < (9.0 * start_x); x += delta_x) {
      x *= (-1.0);
      EXPECT_EQ(1, S21EqDouble(x, x + shifting_koef * local_eps));
      EXPECT_EQ(0, S21EqDouble(x, x + local_eps));

      x *= (-1.0);
      EXPECT_EQ(1, S21EqDouble(x, x + shifting_koef * local_eps));
      EXPECT_EQ(0, S21EqDouble(x, x + local_eps));
    }
    start_x *= 10.0;
  }
  EXPECT_EQ(1, 1);
}

TEST(IsEqMatrixTest, TestEqSize) {
  S21Matrix A(2, 2);
  S21Matrix B(2, 2);

  EXPECT_TRUE(A == B);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  B = A;

  EXPECT_TRUE(B == A);

  B(0, 0) = 1e23;

  EXPECT_FALSE(A == B);

  S21Matrix C;

  EXPECT_FALSE(A == C);
}

TEST(MatrixSumTest, SimpleTest) {
  S21Matrix A(2, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  S21Matrix s21_res = A + A;

  S21Matrix res(2, 2);

  res(0, 0) = 0.0;
  res(0, 1) = 0.2;
  res(1, 0) = 2.0;
  res(1, 1) = 2.2;

  EXPECT_TRUE(res == s21_res);
  EXPECT_FALSE(res == A);

  S21Matrix not_right_res(2, 2);
  EXPECT_FALSE(not_right_res == s21_res);
}

TEST(MatrixSumTest, SumSelfTest) {
  S21Matrix A(2, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  S21Matrix s21_res = A + A;

  A += A;

  EXPECT_TRUE(A == s21_res);

  S21Matrix not_right_res(2, 2);
  EXPECT_FALSE(not_right_res == s21_res);
}

TEST(MatrixSubTest, SubSelfTest) {
  S21Matrix A(2, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  S21Matrix s21_res = A - A;

  S21Matrix res(2, 2);

  EXPECT_TRUE(res == s21_res);
  EXPECT_FALSE(res == A);

  S21Matrix not_right_res(2, 2);
  not_right_res(0, 0) = 11.1;
  EXPECT_FALSE(not_right_res == s21_res);
}

TEST(MatrixSubTest, SimpleTest) {
  S21Matrix A(2, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  S21Matrix res = A;
  S21Matrix s21_res = A + A;
  s21_res -= A;

  EXPECT_TRUE(res == s21_res);
}

TEST(MatrixMulTest, SimpleTest) {
  S21Matrix A(2, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  S21Matrix s21_res = A * A;

  S21Matrix res(2, 2);

  res(0, 0) = 0.10;
  res(0, 1) = 0.11;
  res(1, 0) = 1.10;
  res(1, 1) = 1.31;

  EXPECT_TRUE(res == s21_res);
}

TEST(MatrixMulTest, MulSelfTest) {
  S21Matrix A(2, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  S21Matrix s21_res = A;
  s21_res *= A;

  S21Matrix res(2, 2);

  res(0, 0) = 0.10;
  res(0, 1) = 0.11;
  res(1, 0) = 1.10;
  res(1, 1) = 1.31;

  EXPECT_TRUE(res == s21_res);
}

TEST(MatrixMulNumberTest, SimpleTest) {
  S21Matrix A(2, 2);
  for (int i = 0; i < 4; ++i) {
    A((i / 2), (i % 2)) = i;
  }
  double num = 1.0;

  S21Matrix right_res = A * num;
  S21Matrix left_res = num * A;

  EXPECT_TRUE(right_res == left_res);
  EXPECT_TRUE(right_res == A);
  EXPECT_TRUE(A == left_res);
}

TEST(TransposeTest, SquareTest) {
  S21Matrix A;
  A(0, 0) = 1.0;
  S21Matrix res = A.Transpose();

  EXPECT_TRUE(res == A);

  A.S21ResizeMatrix(2, 2);
  res.S21ResizeMatrix(2, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  res(0, 0) = 0.0;
  res(0, 1) = 1.0;
  res(1, 0) = 0.1;
  res(1, 1) = 1.1;

  A = A.Transpose();

  EXPECT_TRUE(res == A);
}

TEST(TransposeTest, NotSquareTest) {
  S21Matrix A(1, 4);
  S21Matrix res(4, 1);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(0, 2) = 0.2;
  A(0, 3) = 0.3;

  res(0, 0) = 0.0;
  res(1, 0) = 0.1;
  res(2, 0) = 0.2;
  res(3, 0) = 0.3;

  A = A.Transpose();
  EXPECT_TRUE(res == A);

  A.S21ResizeMatrix(2, 3);
  res.S21ResizeMatrix(3, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(0, 2) = 0.2;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;
  A(1, 2) = 1.2;

  res(0, 0) = 0.0;
  res(0, 1) = 1.0;
  res(1, 0) = 0.1;
  res(1, 1) = 1.1;
  res(2, 0) = 0.2;
  res(2, 1) = 1.2;

  A = A.Transpose();

  EXPECT_TRUE(res == A);
}

TEST(Determinant_test, SimpleTest) {
  S21Matrix A;
  double s21_det = A.Determinant();
  double det = 0.0;

  EXPECT_TRUE(S21EqDouble(det, s21_det));

  A.S21ResizeMatrix(2, 2);

  A(0, 0) = 0.0;
  A(0, 1) = 0.1;
  A(1, 0) = 1.0;
  A(1, 1) = 1.1;

  s21_det = A.Determinant();
  det = -0.1;

  EXPECT_TRUE(S21EqDouble(det, s21_det));

  A.S21ResizeMatrix(3, 3);

  A(0, 0) = 0.1;
  A(0, 1) = 0.0;
  A(0, 2) = 0.0;
  A(1, 0) = 0.0;
  A(1, 1) = 1.0;
  A(1, 2) = 0.0;
  A(2, 0) = 0.0;
  A(2, 1) = 0.0;
  A(2, 2) = 1.1;

  s21_det = A.Determinant();
  det = 0.11;

  EXPECT_TRUE(S21EqDouble(det, s21_det));

  A.S21ResizeMatrix(4, 4);

  A(0, 0) = 1.0;
  A(0, 1) = 0.0;
  A(0, 2) = 0.0;
  A(0, 3) = 0.0;
  A(1, 0) = 0.0;
  A(1, 1) = 1.0;
  A(1, 2) = 0.0;
  A(1, 3) = 0.0;
  A(2, 0) = 0.0;
  A(2, 1) = 0.0;
  A(2, 2) = 1.0;
  A(2, 3) = 0.0;
  A(3, 0) = 0.0;
  A(3, 1) = 0.0;
  A(3, 2) = 0.0;
  A(3, 3) = 1.0;

  s21_det = A.Determinant();
  det = 1.0;

  EXPECT_TRUE(S21EqDouble(det, s21_det));
}

TEST(InverseMatrix_test, SimpleTest) {
  S21Matrix A(3, 3);
  A(0, 0) = 1.0;
  A(1, 1) = 1.0;
  A(2, 2) = 1.0;
  S21Matrix B = A.InverseMatrix();

  EXPECT_TRUE(A == B);

  A(0, 0) = 1.0;
  A(1, 1) = 2.0;
  A(2, 2) = 3.0;

  S21Matrix s21_result(3, 3);

  s21_result(0, 0) = 1.0;
  s21_result(1, 1) = 0.5;
  s21_result(2, 2) = 1.0 / 3.0;

  B = A.InverseMatrix();
  EXPECT_TRUE(B == s21_result);

  A.S21ResizeMatrix(1, 1);
  A(0, 0) = 10.0;

  B = A.InverseMatrix();
  s21_result.S21ResizeMatrix(1, 1);
  s21_result(0, 0) = 0.1;

  EXPECT_TRUE(s21_result == B);
}

TEST(Getter_test, SimpleTest) {
  S21Matrix A(1, 1);

  EXPECT_EQ(1, A.GetRows());
  EXPECT_EQ(1, A.GetCols());

  A.S21ResizeMatrix(2, 3);

  EXPECT_EQ(2, A.GetRows());
  EXPECT_EQ(3, A.GetCols());

  A.S21ResizeMatrix(3, 2);

  EXPECT_EQ(3, A.GetRows());
  EXPECT_EQ(2, A.GetCols());
}

TEST(SetRows_test, SimpleTest) {
  S21Matrix A(3, 3);
  for (int i = 0; i < A.GetRows(); ++i) {
    for (int j = 0; j < A.GetCols(); ++j) {
      A(i, j) = (double)(i + j);
    }
  }

  S21Matrix result(2, 3);

  result(0, 0) = 0.0;
  result(0, 1) = 1.0;
  result(0, 2) = 2.0;
  result(1, 0) = 1.0;
  result(1, 1) = 2.0;
  result(1, 2) = 3.0;

  A.SetRows(2);
  EXPECT_TRUE(result == A);
}

TEST(SetCols_test, SimpleTest) {
  S21Matrix A(3, 3);
  for (int i = 0; i < A.GetRows(); ++i) {
    for (int j = 0; j < A.GetCols(); ++j) {
      A(i, j) = (double)(i + j);
    }
  }

  S21Matrix result(3, 2);

  result(0, 0) = 0.0;
  result(0, 1) = 1.0;
  result(1, 0) = 1.0;
  result(1, 1) = 2.0;
  result(2, 0) = 2.0;
  result(2, 1) = 3.0;

  A.SetCols(2);
  S21Matrix s21_result(std::move(A));
  EXPECT_TRUE(result == s21_result);
}

TEST(SetColsRows_test, SimpleTest) {
  S21Matrix A(3, 3);
  for (int i = 0; i < A.GetRows(); ++i) {
    for (int j = 0; j < A.GetCols(); ++j) {
      A(i, j) = (double)(i + j);
    }
  }

  S21Matrix result(2, 2);

  result(0, 0) = 0.0;
  result(0, 1) = 1.0;
  result(1, 0) = 1.0;
  result(1, 1) = 2.0;

  A.SetRows(2);
  A.SetCols(2);
  EXPECT_TRUE(result == A);
}

TEST(Constructors_test, ExeptionsTest) {
  EXPECT_THROW(S21Matrix A(-1, 1), std::out_of_range);
  EXPECT_THROW(S21Matrix A(1, -1), std::out_of_range);
  EXPECT_THROW(S21Matrix A(0, 0), std::out_of_range);
}

TEST(SetRowsCols_test, ExeprionsTest) {
  S21Matrix A;
  EXPECT_THROW(A.SetRows(-1), std::out_of_range);
  EXPECT_THROW(A.SetCols(-1), std::out_of_range);
  EXPECT_THROW(A.S21ResizeMatrix(0, 1), std::out_of_range);
}

TEST(Operators_test, ExeptionsTest) {
  S21Matrix A(2, 2);
  S21Matrix B(3, 3);
  S21Matrix C(4, 5);

  EXPECT_THROW(A + B, std::logic_error);
  EXPECT_THROW(A - B, std::logic_error);
  EXPECT_THROW(A * B, std::logic_error);

  EXPECT_THROW(C.CalcComplements(), std::logic_error);
  EXPECT_THROW(C.Determinant(), std::logic_error);

  EXPECT_THROW(A.InverseMatrix(), std::runtime_error);
  EXPECT_THROW(A(0, -1), std::out_of_range);
  EXPECT_THROW(A(12, 1), std::out_of_range);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
