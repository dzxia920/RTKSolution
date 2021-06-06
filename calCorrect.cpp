#include <cmath>

#include "calCorrect.h"

//Hopfield模型使用的小函数和常量
//定义标准气象元素，单位依次为℃，mbar，1
const double t0 = 20.0;
const double P0 = 1013.25;
const double RH0 = 0.5;
double ComputeEs(double RH, double Ts);
void ComputeMeteoElem(double& t, double& P, double& RH, double h);
double HopfieldCals(double E, double Ps, double es, double Ts, double hs);


void correctEarthRot(sat_t & sat)
{
	double delta_x = OMGE * sat.ttr * sat.pos[1];
	double delta_y = -OMGE * sat.ttr * sat.pos[0];
	sat.pos[0] += delta_x;
	sat.pos[1] += delta_y;
}


double ComputeEs(double RH, double Ts)
{
	double es;
	es = RH*exp(-37.2465 + 0.213166*Ts - 0.000256908*Ts*Ts);
	return es;
}



void ComputeMeteoElem(double& t, double& P, double& RH, double h)
{
	t = t0 - 0.0065*h;
	P = P0*pow((1 - 0.0000266*h), 5.225);
	RH = RH0*exp(-0.0006396*h);
}


double HopfieldCals(double E, double Ps, double es, double Ts, double hs)
{
	//计算辅助量
	double hw = 11000.0;
	double hd = 40136 + 148.72*(Ts - 273.16);
	double Kw = 155.2 * 4810 * es*(hw - hs)*1e-7 / (Ts*Ts);
	double Kd = 155.2*1e-7*Ps*(hd - hs) / Ts;
	double delta_sd = Kd / sin(sqrt(DEG2RAD*E*DEG2RAD*E + DEG2RAD*2.5*DEG2RAD*2.5));
	double delta_sw = Kw / sin(sqrt(DEG2RAD*E*DEG2RAD*E + DEG2RAD*1.5*DEG2RAD*1.5));

	double delta_s = delta_sd + delta_sw;
	return delta_s;
}

double Hopfield(double E, double hs)
{
	double t = 0.0;
	double P = 0.0;
	double RH = 0.0;
	ComputeMeteoElem(t, P, RH, hs);
	double Ts = t + 273.15;
	double es = ComputeEs(RH, Ts);
	double delta_s = HopfieldCals(E, P, es, Ts, hs);
	return delta_s;
}
