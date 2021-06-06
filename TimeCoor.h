#pragma once

/*
	���ļ����������ʱ�䡢����ϵͳ�Ķ��弰ת��

*/

class CommonTime;
class GPST;

//�������� ����������+������sod���
class MJD
{
public:

	//����������ͨ��ʱת��
	CommonTime ToCommonTime();	// ===To Complete 
	GPST ToGPST();				// ===To Complete
	int Days;	//mjd����������
	double sod;	//������
};

class CommonTime
{
public:
	CommonTime() = default;										//Ĭ�Ϲ��캯��
	CommonTime(int y,int m,int d,int h,int min,double s):
		Year(y),Month(m),Day(d),Hour(h),Minute(min),Second(s) {}
	~CommonTime() {}
	MJD ToMJD();												//ͨ��ʱ���������ת��,CHECKED!


	int Year;
	int Month;
	int Day;
	int Hour;
	int Minute;
	double Second;
};

class GPST
{
public:
	GPST() = default;
	~GPST(){}
	
	int week;
	double sow;
};