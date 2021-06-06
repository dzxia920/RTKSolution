#pragma once
#include "global.h"
#include "rinex.h"

/*
*	各项改正说明：
*	由于计算信号发射时卫星位置以及改正了相对论效应，故此处不再改
*	对流层延迟：Hopfield模型
*	电离层延迟：采用消电离层的双频伪距组合
*	地球自转改正：采用P129 相当于卫星坐标旋转了一个角度alpha，
*	TGD：还没改
*/

void correctEarthRot(sat_t& sat);


//给定高度角°，测站高程，利用Hopfield模型计算对流层延迟改正
double Hopfield(double E, double hs);
