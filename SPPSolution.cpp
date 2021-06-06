#include "SPPSolution.h"

const int MAXSOLUTIONSAT = 50;

double VectDot(int m, int n, const double A[], const double B[]);
double Norm(const int n, const double A[]);

double CalDistance(double x[3], double y[3])
{
	double D2 = 0.0;
	for (int i = 0; i < 3; i++)
	{
		D2 += SQR(x[i] - y[i]);
	}
	return sqrt(D2);
}

//���㷽λ�� *************************To be Complete!!!!**********************
double CalAzimuth(double satpos[3],double rcvpos[3])
{
	double azi = 0.0;
	return 0.0;
}

//����߶Ƚ�
double CalElevation(double SatPos[3], double RCVPos[3])
{
	int    i;
	double RcvR, SatRcvR;
	double dPos[3], Elev;       /* ����λ������ջ���λ�ò�ֵ */

	for (i = 0; i<3; i++)
	{
		dPos[i] = SatPos[i] - RCVPos[i];
	}

	RcvR = Norm(3, RCVPos);
	SatRcvR = Norm(3, dPos);

	if (fabs(RcvR * SatRcvR) < 1.0) /* ��������λ�û���ջ�λ���Ƿ�Ϊ0 */
	{
		Elev = PI / 2.0;
	}
	else
	{
		Elev = VectDot(3, 3, RCVPos, dPos) / (RcvR * SatRcvR);
		Elev = PI / 2.0 - acos(Elev);
	}

	return (Elev);
}

void LeastSquare(obs_t & obst, sat_t sat[MAXSATNUM], sol_t & sol)
{
	Matrix B(MAXSOLUTIONSAT, 5);
	Matrix l(MAXSOLUTIONSAT, 5);
	Matrix P = Matrix::Identity(MAXSOLUTIONSAT);
	Matrix x(5, 1);
	Matrix o(MAXSOLUTIONSAT, 1);
	Matrix c(MAXSOLUTIONSAT, 1);
	//Matrix<double, MAXSOLUTIONSAT, 5> B = Matrix<double, MAXSOLUTIONSAT, 5>::zeros();
	//Matrix<double, MAXSOLUTIONSAT, 1> l = Matrix<double, MAXSOLUTIONSAT, 1>::zeros();
	//Matrix<double, MAXSOLUTIONSAT, MAXSOLUTIONSAT> P = Matrix<double, MAXSOLUTIONSAT, MAXSOLUTIONSAT>::Identity();
	//Matrix<double, 5, 1> x = Matrix<double, 5, 1>::zeros();

	//Matrix<double, MAXSOLUTIONSAT, 1> o = Matrix<double, MAXSOLUTIONSAT, 1>::zeros();
	//Matrix<double, MAXSOLUTIONSAT, 1> c = Matrix<double, MAXSOLUTIONSAT, 1>::zeros();

	int numOfObs = 0;
	for (int i = 0; i < MAXSATNUM; i++)
	{		
		if (obst.obsd[i].valid && sat[i].valid) //����ù۲�ֵ��Ч�Ҽ������������λ�ã�����α��
		{
			//��P0α�൥�㶨λ

			//�Ƚ��и�������ļ���
			////To Do ���㷽λ��
			sol.elevation = CalElevation(sat[i].pos, sol.pos);
			//obst.obsd[i].correct.trop = Hopfield();	//��Ҫ����߶ȽǺ͸߳�
			correctEarthRot(sat[i]);				// earthRot ����

			//���Ƿ������ҵļ���
			double rho = CalDistance(sat[i].pos, sol.pos); //���ǺͲ�վ֮��ļ��ξ���
			double li = -(sat[i].pos[0] - sol.pos[0]) / rho;
			double mi = -(sat[i].pos[1] - sol.pos[1]) / rho;
			double ni = -(sat[i].pos[2] - sol.pos[2]) / rho;

			//o(numOfObs,1) = obst.obsd[i].P[0] - obst.obsd[i].correct.trop;
			o(numOfObs, 0) = obst.obsd[i].P[0];
			//o(numOfObs, 0) = 2.54573* obst.obsd[i].P[0] - 1.54573* obst.obsd[i].P[1];
			c(numOfObs, 0) = CalDistance(sat[i].pos,sol.pos); //ע�⣬��ʼʱsol.posΪ000
			B(numOfObs, 0) = li;
			B(numOfObs, 1) = mi;
			B(numOfObs, 2) = ni;
			if (obst.obsd[i].sys == 'G')
			{
				B(numOfObs, 3) = 1;
				B(numOfObs, 4) = 0;
			}
			else if (obst.obsd[i].sys == 'C')
			{
				B(numOfObs, 3) = 0;
				B(numOfObs, 4) = 1;
			}
			numOfObs++;
		}
	}
	sol.validSatNum = numOfObs;
	l = o - c;
	x = (B.transpose()*P*B).inverse()*B.transpose()*P*l;
	sol.pos[0] += x(0, 0);
	sol.pos[1] += x(1, 0);
	sol.pos[2] += x(2, 0);
	sol.dt_GPS = x(3, 0);
	sol.dt_BDS = x(4, 0);
}



/****************************************************************************
Norm

Ŀ�ģ����������ľ���
��ţ�01020

����:
m      A������Ԫ�ظ���
A      ����

����ֵ��Val    ������
****************************************************************************/
double Norm(const int n, const double A[])
{
	if (n <= 0)
	{
		printf("Error dimension in Norm!\n");
		exit(EXIT_FAILURE);
	}

	return (sqrt(VectDot(n, n, A, A)));
}


/****************************************************************************
VectDot    a = A . B

Ŀ�ģ��������������ĵ��
��ţ�01013

����:
m      A������Ԫ�ظ���
n      B������Ԫ�ظ���, Ҫ��m=n
����ֵ��Val    ���
****************************************************************************/
double VectDot(int m, int n, const double A[], const double B[])
{
	int i;
	double Val = 0.0;

	if ((m != n) || (m <= 0) || (n <= 0))
	{
		printf("Error dimension in VectDot!\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i<m; i++)
	{
		Val = Val + *(A + i) * *(B + i);
	}

	return (Val);
}