#pragma once
#include "calCorrect.h"
#include "global.h"
//#include "Matrix.hpp"
#include"MatrixD.h"
#include "rinex.h"

/*
	DRSB 频率间偏差 如果GPS 和BDS 一起
	BDS 和 GPS 看成两个系统，分别选取参考星

	基准站和流动站分别计算卫星位置。
	筛选：截止高度角，如10deg
	可以一步从非差到双差

	参考星选取高度角最大的卫星
	北斗的GEO卫星需要额外处理
*/

//传入该历元的观测值,卫星位置计算结果,解算结果sol
void LeastSquare(obs_t &obst, sat_t sat[MAXSATNUM], sol_t& sol);

//计算高度角
double CalElevation(double SatPos[3], double RCVPos[3]);

//计算两点空间距离
double CalDistance(double x[3], double y[3]);