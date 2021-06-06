#pragma once
// obs and nav �ṹ�嶨��
#include <string>

#include "global.h"
#include "MatrixD.h"
#include "TimeCoor.h"

//������
struct correct_t
{
	double trop;	//���������
	double earthRot;//������ת����
	double iono;	//��������
	double rela;	//�����ЧӦ����
	double TGD;		//Ⱥ�ӳ�

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
	//˫Ƶ��˫ϵͳ
	bool valid;
	char sys;			//ϵͳ��ʶ��G/C
	int prn;			//����PRN��
	double P[2];		//˫Ƶα�����ֵ
	double L[2];		//˫Ƶ�ز���λ����ֵ
	double D[2];		//�����ղ���ֵ
	double SNR[2];		//�����
	double elev;		//�߶Ƚ� ��λ��deg
	//������
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
	MJD mjd;			//ʱ�䣬��MJD��ʾ
	CommonTime ct;
	int n;				//�洢�����ǹ۲�ֵ��
	obsd_t obsd[MAXSATNUM];		//�۲�ֵ�������Ԫ��

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
	std::string obsType[MAXOBSSYSTEM];	//˳����GCERJISM
	double interval;					//���ջ��������
	double approxpos_rcv_xyz[3];		//���ջ�����λ��
	//�ƻ��Ƿ�������˳��
	//GPS�� C1C L1C D1C S1C C2W L2W D2W S2W
	//BDS��	C2I L2I D2I S2I C6I L6I D6I S6I

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



//�����ṹ�壬���ĳ������ĳ����Ԫ������������Ӧ�öԱ����ļ��ԸĽ���
struct eph_t
{
	//���˸���ʱ��ϵͳ��һ���⣬���������ephemeris����һ�µ�
	bool valid;														//��ǰ�����Ƿ���Ч
	char sys;														//ϵͳ��
	int prn;														//PRN��
	double ttr;														//����ʱ��
	MJD toc;														//toc,toe����week+sow,trans time
	int week;														//GPS Week
	double sow;														//������
	double a0, a1, a2;												//���ǵ���ƫ�����١���Ư
	double crc, crs, cuc, cus, cic, cis;							//�����㶯������
	double sqrtA, e, i0, omega0, omega, omegaDot, M0, delta_n, idot;	//�Ų���

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

//�洢���ǵ�λ�ú��ٶ� sat[PRN] ÿ����Ԫ
//�о�Ӧ���ò�վ���洢ÿ�����ǵĴ���ʱ��
struct sat_t
{
	bool valid;			//��Чλ
	char sys;			//G/C
	int prn;			//PRN
	MJD transmitTime;	//��λ�õ�ʱ�䣬��ָ�����źŷ���ʱ��
	double satclk;		//�����Ӳ�
	double vel[3];		//�ٶ�
	double pos[3];		//λ��xyz
	double ttr;			//����ʱ�䣬������תҪ��
	
	sat_t()
	{
		valid = false;
	}
};

//��ŵ���Ԫ������������
struct sol_t
{
	GPST time;		//����ʱ��
	double pos[3];
	double vel[3];
	double dops[4];	//gdop pdop hdop vdop
	int validSatNum; //��Ч������
	double dt_GPS;   //GPS ���ջ��Ӳ�
	double dt_BDS;	//BDS ���ջ��Ӳ�
	double elevation; //�߶Ƚ�
	double azimuth;		//��λ��
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

//���RTK�������Ľṹ��
struct RtkFix
{
	//double pos[3];			//����վ����xyz����
	double _baseline_fix[3];	//���߹̶���
	double ratio;				//�̶����ratio
	Matrix Qaa;					//�̶���ķ���
	//double azim[MAXSATNUM];
	//double elev[MAXSATNUM];
};

//���RTK����������
struct RtkFloat
{
	double pos[3];				//����վ����xyz����
	double _baseline_float[3];	//���߸����
	Matrix x;					//��С���˽������xyz��ģ��������
	Matrix Q;
};




