#pragma once

#include <stdio.h>
#include <stdexcept>

using namespace std;

//最简单的矩阵类
class Matrix
{
public:
	int row;
	int col;
	double **m;

	//构造函数
	Matrix(int rows, int cols);
	//默认构造函数
	Matrix()
	{
		row = 0;
		col = 0;
		m = nullptr;
	}
	//析构函数
	~Matrix();
	//拷贝构造函数
	Matrix(const Matrix& __temp);

	int GetRows()
	{
		return this->row;
	}
	int GetCols()
	{
		return this->col;
	}
	//运算符重载
	//
	//函数调用运算符，可以通过小括号()访问数组里的元素,行列从0开始
	inline double& operator()(int Row, int Col);

	//void zeros()
	//{
	//	for (int i = 0; i < row; i++)
	//		for (int j = 0; j < col; j++)
	//			m[i][j] = 0.0;
	//}

	Matrix& operator=(const Matrix &rhs);

	static Matrix Identity(int rows);


	//运算符重载,其中算术运算符未定义成成员函数
	Matrix& operator+=(const Matrix &rhs);

	Matrix& operator-=(const Matrix &rhs);

	//矩阵相加，矩阵类型必须一致
	Matrix& add(Matrix& M2);

	Matrix transpose();
	Matrix inverse();

	void CopyToArray(double arr[]);

	//成员函数，打印
	void Print();
	void Print2File(FILE *pfile);
};

Matrix operator*(Matrix &lhs, Matrix &rhs);
Matrix operator+(const Matrix& lhs, Matrix& rhs);
Matrix operator-(const Matrix& lhs, Matrix& rhs);

