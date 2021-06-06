#pragma once

//#include<iostream>
#include<initializer_list>

using namespace std;
 

template<class T,int RowsAtCompileTime,int ColsAtCompileTime>
class Matrix
{
public:
	//����Ҫ�õ���̬�����ڴ档��Ϊ��Ԥ�������ά��
	T m_aptr [RowsAtCompileTime][ColsAtCompileTime];
	
	//���캯��
	Matrix();
	//��������
	~Matrix()
	{
	}
	

	//ʹ��initializer_list����ʼ��Matrix
	Matrix(initializer_list<T> ilist)
	{
		if (ilist.size() != RowsAtCompileTime*ColsAtCompileTime)
			throw runtime_error("Matrix Assignment parameters not Accord");
		else {
			auto iter = ilist.begin();
			for (int i = 0; i < RowsAtCompileTime; i++)
			{
				for (int j = 0; j < ColsAtCompileTime; j++)
				{
					m_aptr[i][j] = *iter++;
				}
			}
		}
	}
	//ʹ�������ʼ��
	Matrix(T arr[RowsAtCompileTime][ColsAtCompileTime])
	{
		
		for (int i = 0; i < RowsAtComplieTime; i++)
		{
			for (int j = 0; j < ColsAtComplieTime; j++)
			{
				m_aptr[i][j] = arr[i][j];
			}
		}
	}
	//�������캯��
	Matrix(const Matrix& origin)
	{
		for (int i = 0; i < RowsAtCompileTime; i++)
		{
			for (int j = 0; j < ColsAtCompileTime; j++)
			{
				m_aptr[i][j] = origin.m_aptr[i][j];
			}
		}
	}


	//���������
	//�������
	Matrix<T, RowsAtCompileTime, ColsAtCompileTime>& add(Matrix<T, RowsAtCompileTime, ColsAtCompileTime>& M2);
	
	//��������
	Matrix<T, RowsAtCompileTime, ColsAtCompileTime> inverse();
	//����ת��
	Matrix<T, ColsAtCompileTime, RowsAtCompileTime> transpose();
	
	
	//���������,�������������δ����ɳ�Ա����
	Matrix& operator+=(const Matrix &rhs)
	{
		for (int i = 0; i < RowsAtCompileTime; i++)
		{
			for (int j = 0; j < ColsAtCompileTime; j++)
			{
				m_aptr[i][j] += rhs.m_aptr[i][j];
			}
		}
		return *this;
	}
	Matrix& operator-=(const Matrix &rhs)
	{
		for (int i = 0; i < RowsAtCompileTime; i++)
		{
			for (int j = 0; j < ColsAtCompileTime; j++)
			{
				m_aptr[i][j] -= rhs.m_aptr[i][j];
			}
		}
		return *this;
	}
	
	//������ֵ�����
	Matrix& operator=(const Matrix &rhs)
	{
		for (int i = 0; i < RowsAtCompileTime; i++)
		{
			for (int j = 0; j < ColsAtCompileTime; j++)
			{
				m_aptr[i][j] = rhs.m_aptr[i][j];
			}
		}
		return *this;
	}
	//�������������������ͨ��С����()�����������Ԫ��,���д�0��ʼ
	T& operator()(int Row, int Col)
	{
		return this->m_aptr[Row][Col];
	}
	//�õ�ȫ0����
	static Matrix<T, RowsAtCompileTime, ColsAtCompileTime> zeros()
	{
		Matrix<T, RowsAtCompileTime, ColsAtCompileTime> zeros;
		for (int i = 0; i < RowsAtCompileTime; i++)
			for (int j = 0; j < ColsAtCompileTime; j++)
				zeros(i, j) = 0.0;
		return zeros;
	}
	//�õ���λ��
	static Matrix<T, RowsAtCompileTime, ColsAtCompileTime> Identity()
	{
		Matrix<T, RowsAtCompileTime, ColsAtCompileTime> Identity = Matrix<T, RowsAtCompileTime, ColsAtCompileTime>::zeros();
		for (int i = 0; i < RowsAtCompileTime; i++)
			for (int j = 0; j < ColsAtCompileTime; j++)
				if (i == j)
					Identity(i, j) = 1.0;
		return Identity;
	}
	//��Ա��������ӡ
	void Print();
	
};

//�������
//const�汾��ô���壿��
template<class T, int Rows1, int Rows2,int Cols1,int Cols2>
Matrix<T, Rows1, Cols2> operator*(Matrix<T, Rows1, Cols1>&lhs, Matrix<T, Rows2, Cols2>&rhs)
{
	if (Rows2 != Cols1)
	{
		throw runtime_error("Matrix Multiply must be Row2 equal to Col1");
	}
	Matrix<T, Rows1, Cols2> MatM;
	for (int i = 0; i < Rows1; i++)
	{
		for (int j = 0; j < Cols2; j++)
		{
			MatM(i, j) = 0;
			for (int k = 0; k < Rows2; k++)
			{
				MatM(i, j) += (lhs(i, k)*rhs(k, j));
			}
		}
	}
	return MatM;
}


//�Ӻ����������
template<class C, int Rows, int Cols>
Matrix<C, Rows, Cols> operator+(const Matrix<C, Rows, Cols>&lhs, Matrix<C, Rows, Cols>& rhs)
{
	Matrix<C, Rows, Cols> sum = lhs;
	sum += rhs;
	return sum;
}
//�������������
template<class C, int Rows, int Cols>
Matrix<C, Rows, Cols> operator-(const Matrix<C, Rows, Cols>&lhs, Matrix<C, Rows, Cols>& rhs)
{
	Matrix<C, Rows, Cols> sum = lhs;
	sum -= rhs;
	return sum;
}


//Ĭ�Ϲ��캯����ʵ��
template<class T, int RowsAtCompileTime, int ColsAtCompileTime>
Matrix<T, RowsAtCompileTime, ColsAtCompileTime>::Matrix()
{
	for (int i = 0; i < RowsAtCompileTime; i++)
	{
		for (int j = 0; j < ColsAtCompileTime; j++)
		{
			m_aptr[i][j] = 0;
		}
	}
}

//������ӣ��������ͱ���һ��
template<class T, int RowsAtCompileTime, int ColsAtCompileTime>
Matrix<T, RowsAtCompileTime, ColsAtCompileTime>&
Matrix<T, RowsAtCompileTime, ColsAtCompileTime>::add(Matrix<T, RowsAtCompileTime, ColsAtCompileTime>& M2)
{
	for (int i = 0; i < RowsAtCompileTime; i++)
	{
		for (int j = 0; j < ColsAtCompileTime; j++)
		{
			m_aptr[i][j] += M2.m_aptr[i][j];
		}
	}
	return *this;
}

//�������棬LU�ֽ�
template<class T, int RowsAtCompileTime, int ColsAtCompileTime>
inline Matrix<T, RowsAtCompileTime, ColsAtCompileTime> Matrix<T, RowsAtCompileTime, ColsAtCompileTime>::inverse()
{
	if (RowsAtCompileTime != ColsAtCompileTime)
	{
		throw runtime_error("Matrix Inverse must be rows equal to cols!");
	}
	//n���Ƿ����ά��
	const int n = RowsAtCompileTime;
	Matrix<T, n, n> MatL, MatU, MatInv, MatLInv, MatUInv;
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
				MatU(i,j)=0;
				MatUInv(i, j) = 0;
			}
		}
	}

	//���L��U----------------------------------------
	//���U1j��Li1
	for (int j = 0; j < n; j++)
	{
		MatU(0, j) = this->m_aptr[0][j];
	}
	for (int i = 1; i < n; i++)
	{
		MatL(i, 0) = this->m_aptr[i][0] / MatU(0, 0);
	}
	//��LU����Ԫ��
	for (int r = 1; r < n; r++)
	{
		for (int j = r; j < n; j++)
		{
			MatU(r, j) = this->m_aptr[r][j];
			for (int k = 0; k < r; k++)
			{
				MatU(r, j) -= (MatL(r, k)*MatU(k, j));
			}
		}

		for (int i = r + 1; i < n; i++)
		{
			MatL(i, r) = this->m_aptr[i][r] / MatU(r, r);
			for (int k = 0 ; k < r; k++)
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
	Matrix<T, n, n> MatUTrans = MatU.transpose();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			MatUInv(i, j) = 0;
			
			for (int k = 0; k < i; k++)
			{
				MatUInv(i, j) -= (MatUTrans(i, k)*MatUInv(k, j) / MatUTrans(i, i));
			}
			
		}
		//��ʱ��MatUInv��U��ת�õ���
		MatUInv(i, i) = 1 / MatU(i, i);
	}
	MatUInv = MatUInv.transpose();
	//A�� = U�� * L��
	MatInv = MatUInv*MatLInv;
	//return MatU;
	return MatInv;
}

//����ת��
template<class T, int RowsAtCompileTime, int ColsAtCompileTime>
inline Matrix<T, ColsAtCompileTime, RowsAtCompileTime> Matrix<T, RowsAtCompileTime, ColsAtCompileTime>::transpose()
{
	Matrix<T, ColsAtCompileTime,RowsAtCompileTime> MatTrans;
	for (int i = 0; i < RowsAtCompileTime; i++)
	{
		for (int j = 0; j < ColsAtCompileTime; j++)
		{
			MatTrans(j, i) = this->m_aptr[i][j];
		}
	}
	return MatTrans;
}



//����� ��ӡ����
template<class T, int RowsAtComplieTime, int ColsAtComplieTime>
inline void Matrix<T, RowsAtComplieTime, ColsAtComplieTime>::Print()
{
	for (int i = 0; i < RowsAtComplieTime; i++)
	{
		for (int j = 0; j < ColsAtComplieTime; j++)
		{
			printf("%5f ", this->m_aptr[i][j]);
			//cout << this->m_aptr[i][j] << " ";
		}
		//cout << endl;
		printf("\n");
	}
}