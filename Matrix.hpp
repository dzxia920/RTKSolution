#pragma once

//#include<iostream>
#include<initializer_list>

using namespace std;
 

template<class T,int RowsAtCompileTime,int ColsAtCompileTime>
class Matrix
{
public:
	//不需要用到动态分配内存。因为可预设数组的维度
	T m_aptr [RowsAtCompileTime][ColsAtCompileTime];
	
	//构造函数
	Matrix();
	//析构函数
	~Matrix()
	{
	}
	

	//使用initializer_list来初始化Matrix
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
	//使用数组初始化
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
	//拷贝构造函数
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


	//矩阵的运算
	//矩阵相加
	Matrix<T, RowsAtCompileTime, ColsAtCompileTime>& add(Matrix<T, RowsAtCompileTime, ColsAtCompileTime>& M2);
	
	//矩阵求逆
	Matrix<T, RowsAtCompileTime, ColsAtCompileTime> inverse();
	//矩阵转置
	Matrix<T, ColsAtCompileTime, RowsAtCompileTime> transpose();
	
	
	//运算符重载,其中算术运算符未定义成成员函数
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
	
	//拷贝赋值运算符
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
	//函数调用运算符，可以通过小括号()访问数组里的元素,行列从0开始
	T& operator()(int Row, int Col)
	{
		return this->m_aptr[Row][Col];
	}
	//得到全0矩阵
	static Matrix<T, RowsAtCompileTime, ColsAtCompileTime> zeros()
	{
		Matrix<T, RowsAtCompileTime, ColsAtCompileTime> zeros;
		for (int i = 0; i < RowsAtCompileTime; i++)
			for (int j = 0; j < ColsAtCompileTime; j++)
				zeros(i, j) = 0.0;
		return zeros;
	}
	//得到单位阵
	static Matrix<T, RowsAtCompileTime, ColsAtCompileTime> Identity()
	{
		Matrix<T, RowsAtCompileTime, ColsAtCompileTime> Identity = Matrix<T, RowsAtCompileTime, ColsAtCompileTime>::zeros();
		for (int i = 0; i < RowsAtCompileTime; i++)
			for (int j = 0; j < ColsAtCompileTime; j++)
				if (i == j)
					Identity(i, j) = 1.0;
		return Identity;
	}
	//成员函数，打印
	void Print();
	
};

//矩阵相乘
//const版本怎么定义？？
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


//加号运算符重载
template<class C, int Rows, int Cols>
Matrix<C, Rows, Cols> operator+(const Matrix<C, Rows, Cols>&lhs, Matrix<C, Rows, Cols>& rhs)
{
	Matrix<C, Rows, Cols> sum = lhs;
	sum += rhs;
	return sum;
}
//减号运算符重载
template<class C, int Rows, int Cols>
Matrix<C, Rows, Cols> operator-(const Matrix<C, Rows, Cols>&lhs, Matrix<C, Rows, Cols>& rhs)
{
	Matrix<C, Rows, Cols> sum = lhs;
	sum -= rhs;
	return sum;
}


//默认构造函数的实现
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

//矩阵相加，矩阵类型必须一致
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

//矩阵求逆，LU分解
template<class T, int RowsAtCompileTime, int ColsAtCompileTime>
inline Matrix<T, RowsAtCompileTime, ColsAtCompileTime> Matrix<T, RowsAtCompileTime, ColsAtCompileTime>::inverse()
{
	if (RowsAtCompileTime != ColsAtCompileTime)
	{
		throw runtime_error("Matrix Inverse must be rows equal to cols!");
	}
	//n就是方阵的维数
	const int n = RowsAtCompileTime;
	Matrix<T, n, n> MatL, MatU, MatInv, MatLInv, MatUInv;
	//构建上三角阵U和下三角阵L 以及他们的逆矩阵LInv,UInv
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

	//求解L和U----------------------------------------
	//求解U1j和Li1
	for (int j = 0; j < n; j++)
	{
		MatU(0, j) = this->m_aptr[0][j];
	}
	for (int i = 1; i < n; i++)
	{
		MatL(i, 0) = this->m_aptr[i][0] / MatU(0, 0);
	}
	//求LU后续元素
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

	//LU求逆
	//下三角矩阵L求逆
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
	//上三角矩阵U求逆，转置求逆再转置  !!是对U的转置求逆
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
		//此时的MatUInv是U的转置的逆
		MatUInv(i, i) = 1 / MatU(i, i);
	}
	MatUInv = MatUInv.transpose();
	//A逆 = U逆 * L逆
	MatInv = MatUInv*MatLInv;
	//return MatU;
	return MatInv;
}

//矩阵转置
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



//输出， 打印操作
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