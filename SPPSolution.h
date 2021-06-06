#pragma once
#include "calCorrect.h"
#include "global.h"
//#include "Matrix.hpp"
#include"MatrixD.h"
#include "rinex.h"

/*
	DRSB Ƶ�ʼ�ƫ�� ���GPS ��BDS һ��
	BDS �� GPS ��������ϵͳ���ֱ�ѡȡ�ο���

	��׼վ������վ�ֱ��������λ�á�
	ɸѡ����ֹ�߶Ƚǣ���10deg
	����һ���ӷǲ˫��

	�ο���ѡȡ�߶Ƚ���������
	������GEO������Ҫ���⴦��
*/

//�������Ԫ�Ĺ۲�ֵ,����λ�ü�����,������sol
void LeastSquare(obs_t &obst, sat_t sat[MAXSATNUM], sol_t& sol);

//����߶Ƚ�
double CalElevation(double SatPos[3], double RCVPos[3]);

//��������ռ����
double CalDistance(double x[3], double y[3]);