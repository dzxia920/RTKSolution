#include "MatrixD.h"

inline Matrix::Matrix(int rows, int cols)
{
	row = rows;
	col = cols;
	m = new double*[row];
	for (int i = 0; i < row; i++)
	{
		m[i] = new double[col];
	}

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			m[i][j] = 0.0;
		}
	}
}

inline Matrix::~Matrix()
{
	if (m)
	{
		for (int i = 0; i < row; i++)
		{
			delete[] m[i];
			//m[i] = nullptr;
		}
		delete[] m;
		//m = nullptr;
	}
}

inline Matrix::Matrix(const Matrix& __temp)
{
	//if (m)
	//{
	//	for (int i = 0; i < row; i++)
	//	{
	//		delete[] m[i];
	//		//m[i] = nullptr;
	//	}
	//	delete[] m;
	//	//m = nullptr;
	//}
	row = __temp.row;
	col = __temp.col;
	if (row>0 && col>0)
	{
		m = new double*[row];
		for (int i = 0; i<row; ++i)
			m[i] = new double[col];
		for (int i = 0; i<row; ++i)
			for (int j = 0; j<col; ++j)
				m[i][j] = __temp.m[i][j];
	}
	else
	{
		row = 0;
		col = 0;
		m = NULL;
	}
	return;
}

double& Matrix::operator()(int Row, int Col)
{
	return this->m[Row][Col];
}

Matrix& Matrix::operator=(const Matrix &rhs)
{
	//������
	if (m)
	{
		for (int i = 0; i < row; i++)
		{
			delete[] m[i];
			m[i] = nullptr;
		}
		delete[] m;
		m = nullptr;
	}

	row = rhs.row;
	col = rhs.col;
	//����
	m = new double *[row];
	for (int i = 0; i < row; i++)
	{
		m[i] = new double[col];
	}

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			m[i][j] = rhs.m[i][j];
		}
	}
	return *this;
}

Matrix Matrix::Identity(int rows)
{
	Matrix mat(rows, rows);
	for (int i = 0; i < rows; i++)
	{
		mat(i, i) = 1.0;
	}
	return mat;
}

Matrix& Matrix::operator+=(const Matrix &rhs)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			m[i][j] += rhs.m[i][j];
		}
	}
	return *this;
}

Matrix& Matrix::operator-=(const Matrix &rhs)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			m[i][j] -= rhs.m[i][j];
		}
	}
	return *this;
}

Matrix& Matrix::add(Matrix& M2)
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			m[i][j] += M2.m[i][j];
		}
	}
	return *this;
}

void Matrix::Print()
{
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			printf("% f ", this->m[i][j]);
		}
		printf("\n");
	}
}

void Matrix::Print2File(FILE *pfile)
{
	
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			fprintf(pfile,"% f ", this->m[i][j]);
		}
		fprintf(pfile,"\n");
	}
	fprintf(pfile, "\n");


}

void Matrix::CopyToArray(double arr[])
{

	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			m[i][j] = arr[i*col + j];
		}
	}
}


//����ĺ���
//�������
Matrix operator*(Matrix &lhs, Matrix &rhs)
{

	if (rhs.row != lhs.col)
	{
		throw runtime_error("Matrix Multiply must be Row2 equal to Col1");
	}
	Matrix MatM(lhs.row, rhs.col);
	for (int i = 0; i < lhs.row; i++)
	{
		for (int j = 0; j < rhs.col; j++)
		{
			MatM(i, j) = 0;
			for (int k = 0; k <rhs.row; k++)
			{
				MatM(i, j) += (lhs(i, k)*rhs(k, j));
			}
		}
	}
	return MatM;
}

//�Ӻ����������
Matrix operator+(const Matrix& lhs, Matrix& rhs)
{
	Matrix sum(lhs.row, lhs.col);
	sum = lhs;
	sum += rhs;
	return sum;
}

//�������������
Matrix operator-(const Matrix& lhs, Matrix& rhs)
{
	Matrix sum(lhs.row, lhs.col);
	sum = lhs;
	sum -= rhs;
	return sum;
}

//����ת��
Matrix Matrix::transpose()
{
	Matrix MatTrans(col, row);
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			MatTrans(j, i) = this->m[i][j];
		}
	}
	return MatTrans;
}

//�������棬LU�ֽ�
Matrix Matrix::inverse()
{
	if (row != col)
	{
		throw runtime_error("Matrix Inverse must be rows equal to cols!");
	}
	//n���Ƿ����ά��
	const int n = row;
	Matrix  MatL(n, n), MatU(n, n), MatInv(n, n), MatLInv(n, n), MatUInv(n, n);
	//������������U����������L �Լ����ǵ������LInv,UInv
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				MatL(i, j) = 1;
				MatLInv(i, j) = 1;
			}
			else if (i < j)
			{
				MatL(i, j) = 0;
				MatLInv(i, j) = 0;
			}
			else {
				MatU(i, j) = 0;
				MatUInv(i, j) = 0;
			}
		}
	}

	//���L��U----------------------------------------
	//���U1j��Li1
	for (int j = 0; j < n; j++)
	{
		MatU(0, j) = this->m[0][j];
	}
	for (int i = 1; i < n; i++)
	{
		MatL(i, 0) = this->m[i][0] / MatU(0, 0);
	}
	//��LU����Ԫ��
	for (int r = 1; r < n; r++)
	{
		for (int j = r; j < n; j++)
		{
			MatU(r, j) = this->m[r][j];
			for (int k = 0; k < r; k++)
			{
				MatU(r, j) -= (MatL(r, k)*MatU(k, j));
			}
		}

		for (int i = r + 1; i < n; i++)
		{
			MatL(i, r) = this->m[i][r] / MatU(r, r);
			for (int k = 0; k < r; k++)
			{
				MatL(i, r) -= MatL(i, k)*MatU(k, r) / MatU(r, r);
			}
		}
	}
	//-----------------------------------------------------------

	//LU����
	//�����Ǿ���L����
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			MatLInv(i, j) = 0;
			for (int k = 0; k < i; k++)
			{
				MatLInv(i, j) -= (MatL(i, k)*MatLInv(k, j) / MatL(i, i));
			}
		}
	}
	//�����Ǿ���U���棬ת��������ת��  !!�Ƕ�U��ת������
	Matrix MatUTrans(n, n);
	MatUTrans = MatU.transpose();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			MatUInv(i, j) = 0;

			for (int k = 0; k < i; k++)
			{
				MatUInv(i, j) -= ((MatUTrans(i, k)*MatUInv(k, j)) / MatUTrans(i, i));
			}

		}
		//��ʱ��MatUInv��U��ת�õ���
		MatUInv(i, i) = 1.0 / MatU(i, i);
	}
	MatUInv = MatUInv.transpose();
	//A�� = U�� * L��
	MatInv = MatUInv*MatLInv;
	//return MatU;
	return MatInv;
}
