#pragma once
// obs and nav 结构体定义
#include <string>

#include "global.h"
#include "MatrixD.h"
#include "TimeCoor.h"

//改正项
struct correct_t
{
	double trop;	//对流层改正
	double earthRot;//地球自转改正
	double iono;	//电离层改正
	double rela;	//相对论效应改正
	double TGD;		//群延迟

	correct_t()
	{
		trop = 0.0;
		earthRot = 0.0;
		iono = 0.0;
		rela = 0.0;
		TGD = 0.0;
	}
};

struct obsd_t
{
	//双频、双系统
	bool valid;
	char sys;			//系统标识符G/C
	int prn;			//卫星PRN号
	double P[2];		//双频伪距测量值
	double L[2];		//双频载波相位测量值
	double D[2];		//多普勒测量值
	double SNR[2];		//信噪比
	double elev;		//高度角 单位：deg
	//改正项
	correct_t correct;

	obsd_t()
	{
		elev = 0.0;
		valid = false;
	}
	void Print()
	{
		if (valid)
			printf("%c%2d%14.3f  %14.3f  %14.3f  %14.3f  %14.3f  %14.3f  %14.3f  %14.3f\n",
				sys, prn, P[0], L[0], D[0], SNR[0], P[1], L[1], D[1], SNR[1]);
	}
};

struct obs_t
{
	MJD mjd;			//时间，用MJD表示
	CommonTime ct;
	int n;				//存储的卫星观测值数
	obsd_t obsd[MAXSATNUM];		//观测值数组的首元素

	void Print()
	{
		//CommonTime ct = mjd.ToCommonTime();
		printf("%d %d %d %d %d %f %d\n", ct.Year,ct.Month,ct.Day,ct.Hour,ct.Minute,ct.Second, n);
		for (int i = 0; i < MAXSATNUM; i++)
			obsd[i].Print();
	}
};



class obs_header
{
public:
	std::string obsType[MAXOBSSYSTEM];	//顺序是GCERJISM
	double interval;					//接收机采样间隔
	double approxpos_rcv_xyz[3];		//接收机概略位置
	//计划是返回两个顺序
	//GPS： C1C L1C D1C S1C C2W L2W D2W S2W
	//BDS：	C2I L2I D2I S2I C6I L6I D6I S6I

	obs_header()
	{
		interval = 0.0;
		approxpos_rcv_xyz[0] = 0.0;
		approxpos_rcv_xyz[1] = 0.0;
		approxpos_rcv_xyz[2] = 0.0;
	}
	
	void calobsTypeSequence(int obsTypeSequence[MAXOBSSYSTEM][8])
	{
		//GPS
		for (int i = 0; i < 8; i++)
		{
			//get index 0 1 2 3...
			obsTypeSequence[0][i] = (obsType[0].find(GPSOBSTYPE[i]) - 7) / 4;
			obsTypeSequence[1][i] = (obsType[1].find(BDSOBSTYPE[i]) - 7) / 4;
		}
	}


};



//星历结构体，存放某颗卫星某个历元的星历（后续应该对北斗的加以改进）
struct eph_t
{
	//除了各自时间系统不一致外，其余参数在ephemeris上是一致的
	bool valid;														//当前星历是否有效
	char sys;														//系统号
	int prn;														//PRN号
	double ttr;														//传播时间
	MJD toc;														//toc,toe就是week+sow,trans time
	int week;														//GPS Week
	double sow;														//周内秒
	double a0, a1, a2;												//卫星的钟偏、钟速、钟漂
	double crc, crs, cuc, cus, cic, cis;							//六个摄动改正项
	double sqrtA, e, i0, omega0, omega, omegaDot, M0, delta_n, idot;	//九参数

	eph_t()
	{
		valid = false;
	}
};

struct nav_header
{

};

struct nav_t
{

};

//存储卫星的位置和速度 sat[PRN] 每个历元
//感觉应该让测站来存储每颗卫星的传播时间
struct sat_t
{
	bool valid;			//有效位
	char sys;			//G/C
	int prn;			//PRN
	MJD transmitTime;	//此位置的时间，是指卫星信号发射时刻
	double satclk;		//卫星钟差
	double vel[3];		//速度
	double pos[3];		//位置xyz
	double ttr;			//传播时间，地球自转要用
	
	sat_t()
	{
		valid = false;
	}
};

//存放单历元解算结果的数组
struct sol_t
{
	GPST time;		//解算时间
	double pos[3];
	double vel[3];
	double dops[4];	//gdop pdop hdop vdop
	int validSatNum; //有效卫星数
	double dt_GPS;   //GPS 接收机钟差
	double dt_BDS;	//BDS 接收机钟差
	double elevation; //高度角
	double azimuth;		//方位角
	sol_t()
	{
		for (int i = 0; i < 3; i++)
		{
			pos[i] = 0.0;
			vel[i] = 0.0;
		}
		elevation = 0.0;
		azimuth = 0.0;
	}
};

//存放RTK解算结果的结构体
struct RtkFix
{
	//double pos[3];			//流动站解算xyz坐标
	double _baseline_fix[3];	//基线固定解
	double ratio;				//固定解的ratio
	Matrix Qaa;					//固定解的方差
	//double azim[MAXSATNUM];
	//double elev[MAXSATNUM];
};

//存放RTK浮点解解算结果
struct RtkFloat
{
	double pos[3];				//流动站解算xyz坐标
	double _baseline_float[3];	//基线浮点解
	Matrix x;					//最小二乘解算出的xyz及模糊度向量
	Matrix Q;
};




