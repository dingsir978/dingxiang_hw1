#include "algebra.h"
#include <stdio.h>
#include <math.h>
#define EPSILON 1e-9 // 定义极小值，用于浮点数比较

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
// 检查行列是否一致
        if (a.rows != b.rows || a.cols != b.cols) {
            printf("Error: Matrix a and b must have the same rows and cols.\n");
            return create_matrix(0, 0);
        }
// 创建结果矩阵
        Matrix result = create_matrix(a.rows, a.cols);
        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < a.cols; j++) {
                result.data[i][j] = a.data[i][j] + b.data[i][j];
            }
        }
        return result;
}

Matrix sub_matrix(Matrix a, Matrix b)
{
 // 检查行列是否一致
 if (a.rows != b.rows || a.cols != b.cols) {
    printf("Error: Matrix a and b must have the same rows and cols.\n");
    return create_matrix(0, 0);
}
// 创建结果矩阵
Matrix result = create_matrix(a.rows, a.cols);
for (int i = 0; i < a.rows; i++) {
    for (int j = 0; j < a.cols; j++) {
        result.data[i][j] = a.data[i][j] - b.data[i][j];
    }
}
return result;
}

Matrix mul_matrix(Matrix a, Matrix b)
{
// 检查行列是否一致
if (a.cols != b.rows) {
    printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
    return create_matrix(0, 0);
}
// 创建结果矩阵，行数为a.rows，列数为b.cols
Matrix result = create_matrix(a.rows, b.cols);
for (int i = 0; i < a.rows; i++) {
    for (int j = 0; j < b.cols; j++) {
        // 计算点积
        double sum = 0.0;
        for (int k = 0; k < a.cols; k++) {
            sum += a.data[i][k] * b.data[k][j];
        }
        result.data[i][j] = sum;
    }
}
return result;
}

Matrix scale_matrix(Matrix a, double k)
{
// 创建结果矩阵并乘以标量k
Matrix result = create_matrix(a.rows, a.cols);
for (int i = 0; i < a.rows; i++) {
    for (int j = 0; j < a.cols; j++) {
        result.data[i][j] = a.data[i][j] * k;
    }
}
return result;
}

Matrix transpose_matrix(Matrix a)
{
      // 转置矩阵的行列互换
      Matrix result = create_matrix(a.cols, a.rows);
      for (int i = 0; i < a.cols; i++) {
          for (int j = 0; j < a.rows; j++) {
              result.data[i][j] = a.data[j][i];
          }
      }
      return result;
  }
  

double det_matrix(Matrix a)
{
 // 检查是否为方阵
 if (a.rows != a.cols) {
    printf("Error: The matrix must be a square matrix.\n");
    return 0.0;
}
int n = a.rows;
double det = 1.0;
double temp[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
// 复制数据到临时数组
for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
        temp[i][j] = a.data[i][j];
    }
}
int sign = 1; // 行列式符号

// 高斯消元转换为上三角矩阵
for (int i = 0; i < n; i++) {
    // 寻找主元
    int max_row = i;
    for (int k = i; k < n; k++) {
        if (fabs(temp[k][i]) > fabs(temp[max_row][i])) {
            max_row = k;
        }
    }
    // 处理全零列的情况
    if (fabs(temp[max_row][i]) < EPSILON) {
        return 0.0;
    }
    // 交换行并更新符号
    if (max_row != i) {
        for (int j = 0; j < n; j++) {
            double tmp = temp[i][j];
            temp[i][j] = temp[max_row][j];
            temp[max_row][j] = tmp;
        }
        sign *= -1;
    }
    // 消元
    for (int k = i + 1; k < n; k++) {
        double factor = temp[k][i] / temp[i][i];
        for (int j = i; j < n; j++) {
            temp[k][j] -= factor * temp[i][j];
        }
    }
}
// 计算对角线乘积
for (int i = 0; i < n; i++) {
    det *= temp[i][i];
}
return det * sign;
}

Matrix inv_matrix(Matrix a)
{
// 检查是否为方阵
      if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    int n = a.rows;
    // 检查行列式是否为0
    double det = det_matrix(a);
    if (fabs(det) < EPSILON) {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    // 构造增广矩阵 [A|I]
    double aug[MAX_MATRIX_SIZE][2 * MAX_MATRIX_SIZE] = {0};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i][j] = a.data[i][j];
        }
        aug[i][n + i] = 1.0;
    }
    // 高斯-约当消元
    for (int i = 0; i < n; i++) {
        // 寻找主元
        int max_row = i;
        for (int k = i; k < n; k++) {
            if (fabs(aug[k][i]) > fabs(aug[max_row][i])) {
                max_row = k;
            }
        }
        // 交换行
        if (max_row != i) {
            for (int j = 0; j < 2 * n; j++) {
                double tmp = aug[i][j];
                aug[i][j] = aug[max_row][j];
                aug[max_row][j] = tmp;
            }
        }
        // 归一化当前行
        double pivot = aug[i][i];
        for (int j = 0; j < 2 * n; j++) {
            aug[i][j] /= pivot;
        }
        // 消去其他行
        for (int k = 0; k < n; k++) {
            if (k != i && fabs(aug[k][i]) > EPSILON) {
                double factor = aug[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
    }
    // 提取逆矩阵部分
    Matrix inv = create_matrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv.data[i][j] = aug[i][j + n];
        }
    }
    return inv;
}

int rank_matrix(Matrix a)
{
    int rank = 0;
    int rows = a.rows;
    int cols = a.cols;
    double mat[MAX_MATRIX_SIZE][MAX_MATRIX_SIZE];
    // 复制数据到临时数组
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            mat[i][j] = a.data[i][j];
        }
    }
    // 高斯消元求行阶梯形
    for (int col = 0, row = 0; col < cols && row < rows; col++) {
        int pivot = row;
        // 寻找最大主元
        for (int i = row; i < rows; i++) {
            if (fabs(mat[i][col]) > fabs(mat[pivot][col])) {
                pivot = i;
            }
        }
        if (fabs(mat[pivot][col]) < EPSILON) continue;
        // 交换行
        if (pivot != row) {
            for (int j = col; j < cols; j++) {
                double tmp = mat[row][j];
                mat[row][j] = mat[pivot][j];
                mat[pivot][j] = tmp;
            }
        }
        // 消去下方元素
        for (int i = row + 1; i < rows; i++) {
            double factor = mat[i][col] / mat[row][col];
            for (int j = col; j < cols; j++) {
                mat[i][j] -= factor * mat[row][j];
            }
        }
        rank++;
        row++;
    }
    return rank;
}

double trace_matrix(Matrix a)
{
 // 检查是否为方阵
       if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0.0;
    }
    double trace = 0.0;
    for (int i = 0; i < a.rows; i++) {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}