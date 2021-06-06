#pragma once

/*
	该文件包含所需的时间、坐标系统的定义及转换

*/

class CommonTime;
class GPST;

//简化儒略日 由整数天数+日内秒sod组成
class MJD
{
public:

	//简化儒略日向通用时转化
	CommonTime ToCommonTime();	// ===To Complete 
	GPST ToGPST();				// ===To Complete
	int Days;	//mjd的整数部分
	double sod;	//日内秒
};

class CommonTime
{
public:
	CommonTime() = default;										//默认构造函数
	CommonTime(int y,int m,int d,int h,int min,double s):
		Year(y),Month(m),Day(d),Hour(h),Minute(min),Second(s) {}
	~CommonTime() {}
	MJD ToMJD();												//通用时向简化儒略日转换,CHECKED!


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