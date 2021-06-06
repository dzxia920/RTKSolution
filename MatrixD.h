#pragma once

#include <stdio.h>
#include <stdexcept>

using namespace std;

//��򵥵ľ�����
class Matrix
{
public:
	int row;
	int col;
	double **m;

	//���캯��
	Matrix(int rows, int cols);
	//Ĭ�Ϲ��캯��
	Matrix()
	{
		row = 0;
		col = 0;
		m = nullptr;
	}
	//��������
	~Matrix();
	//�������캯��
	Matrix(const Matrix& __temp);

	int GetRows()
	{
		return this->row;
	}
	int GetCols()
	{
		return this->col;
	}
	//���������
	//
	//�������������������ͨ��С����()�����������Ԫ��,���д�0��ʼ
	inline double& operator()(int Row, int Col);

	//void zeros()
	//{
	//	for (int i = 0; i < row; i++)
	//		for (int j = 0; j < col; j++)
	//			m[i][j] = 0.0;
	//}

	Matrix& operator=(const Matrix &rhs);

	static Matrix Identity(int rows);


	//���������,�������������δ����ɳ�Ա����
	Matrix& operator+=(const Matrix &rhs);

	Matrix& operator-=(const Matrix &rhs);

	//������ӣ��������ͱ���һ��
	Matrix& add(Matrix& M2);

	Matrix transpose();
	Matrix inverse();

	void CopyToArray(double arr[]);

	//��Ա��������ӡ
	void Print();
	void Print2File(FILE *pfile);
};

Matrix operator*(Matrix &lhs, Matrix &rhs);
Matrix operator+(const Matrix& lhs, Matrix& rhs);
Matrix operator-(const Matrix& lhs, Matrix& rhs);

