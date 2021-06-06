#include "TimeCoor.h"
#include <cmath>


CommonTime MJD::ToCommonTime()
{
	// TODO: 在此处插入 return 语句
	
	CommonTime ct(0,0,0,0,0,0.0);

	return ct;
}

GPST MJD::ToGPST()
{
	GPST gpst;
	gpst.week = (Days - 44244) / 7;
	gpst.sow = (Days - 44244 - gpst.week * 7)*86400.0 + sod;
	return gpst;
}

MJD CommonTime::ToMJD()
{
	double y = Year;
	double m = Month;

	if (Month < 3)
	{
		y = y - 1;
		m = m + 12;
	}
	//double h = Hour + Minute / 60.0 + Second / 3600.0;
	//double JD = floor(365.25*y) + floor(30.6001*(m + 1)) + Day + h / 24.0 + 1720981.5;
	//double mjd = JD - 2400000.5;

	
	MJD _mjd;
	_mjd.Days = static_cast<int>(floor(365.25*y) + floor(30.6001*(m + 1)) + Day + 1720981.5 - 2400000.5);
	_mjd.sod = Hour*3600 + Minute*60 + Second;
	return _mjd;
}
