#include "calSatPV.h"

/************函数内部使用的小函数*******************/
double ComputeE(double Mt, double e);	//迭代计算Et

//这里默认传进来的是某颗卫星的一组有效星历
void calSingleSatPV(MJD imputT, const eph_t& eph, sat_t& sat)
{
	// 卫星位置根据啥时间来算，toe？

	//需要将传进来的mjd转成sow 周内秒
	double t = imputT.ToGPST().sow;
	double toe = eph.sow;
	//平均角速度n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(mu) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//计算观测瞬间 t 时刻卫星的平近点角
	//这里的t为周内秒，和参考时刻不同一周时需要使t-toe的绝对值在302400间
	double Mt = eph.M0 + n*(t - toe);
	//判断Mt是否在0-2pi之间
	while (Mt > 2 * PI)
	{
		Mt = Mt - 2 * PI;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI;
	}

	//计算偏近点角
	double Et = ComputeE(Mt,eph.e);

	//计算真近点角
	double upper = sin(Et)*sqrt(1 -eph.e*eph.e);
	double under = cos(Et) - eph.e;
	//double ft = atan(upper / under);
	double ft = atan2(upper, under);

	//计算升交距角（未经改正的）
	double ut_pie = eph.omega + ft;

	//计算卫星向径
	double rt_pie = eph.sqrtA*eph.sqrtA*(1 - eph.e*cos(Et));

	//计算摄动改正项
	double delta_ut = eph.cuc*cos(2 * ut_pie) + eph.cus*sin(2 * ut_pie);
	double delta_rt = eph.crc*cos(2 * ut_pie) + eph.crs*sin(2 * ut_pie);
	double delta_it = eph.cic*cos(2 * ut_pie) + eph.cis*sin(2 * ut_pie);

	//进行摄动改正
	double ut = ut_pie + delta_ut;
	double rt = rt_pie + delta_rt;
	double it = eph.i0 + eph.idot*(t - toe) + delta_it;

	//计算卫星在轨道平面坐标系中的位置
	double xt = rt*cos(ut);
	double yt = rt*sin(ut);

	//计算观测瞬间升交点经度
	double Lt = eph.omega0 + (eph.omegaDot - OMGE)*t - eph.omegaDot*toe;

	//计算卫星在地固坐标系下的坐标
	sat.pos[0] = xt*cos(Lt) - yt*cos(it)*sin(Lt);
	sat.pos[1] = xt*sin(Lt) + yt*cos(it)*cos(Lt);
	sat.pos[2] = yt*sin(it);

	//计算钟差
	//sat.satclk = eph.a0 + eph.a1*(t - toe) + eph.a2*(t - toe)*(t - toe);
	
}

// 计算deltaTr
double CaldeltaTr(MJD imputT, const eph_t& eph)
{
	//需要将传进来的mjd转成sow 周内秒
	double t = imputT.ToGPST().sow;
	double toe = eph.sow;
	//平均角速度n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(mu) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//计算观测瞬间 t 时刻卫星的平近点角
	//这里的t为周内秒，和参考时刻不同一周时需要使t-toe的绝对值在302400间
	double Mt = eph.M0 + n*(t - toe);
	//判断Mt是否在0-2pi之间
	while (Mt > 2 * PI)
	{
		Mt = Mt - 2 * PI;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI;
	}

	//计算偏近点角
	double Et = ComputeE(Mt, eph.e);

	double deltaTr = -2 * sqrt(eph.sqrtA*eph.sqrtA*mu) * eph.e*sin(Et) / CLIGHT / CLIGHT;
	return deltaTr;
}

// input放的应该是obs.mjd，所以sat存放的是那个历元所有卫星的坐标
void CalBatchSatPV(MJD inputT,  eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM])
{	
	//观测值的toe
	double toe_Input = inputT.ToGPST().sow;

	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		//寻找最近的星历，对某颗卫星而言
		for (int j = 0; j < EPOCHNUM; j++)
		{
			if (eph[i][j].valid)
			{
				//对于GPS来说是3600s
				//如果满足星历toe与传入的toe之间相差3600以内，选择该组星历计算卫星位置
				if (abs(eph[i][j].sow - toe_Input) <= 3600.0)
				{
					calSingleSatPV(inputT, eph[i][j], sat[i]);
					sat[i].valid = true;
				}
			}
		}
	}
}

//迭代计算卫星信号发射时刻的卫星位置,感觉还需要伪距测量值。
void calSatTransPV(obs_t& obst, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM])
{
	//根据接收到的时间来确定发射时间

	//梅开二度，需要将所有标志为置成false
	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		sat[i].valid = false;
	}

	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		if (obst.obsd[i].valid)		//判断他有效
		{
			
			//这里可以用双频消电离层伪距组合 Page107 对于GPS和BDS可能不一样吧
			double P_if = 2.54573* obst.obsd[i].P[0] - 1.54573* obst.obsd[i].P[1];
			double toe_Input = obst.mjd.ToGPST().sow;

			//寻找最近的星历，对某颗卫星而言
			for (int j = 0; j < EPOCHNUM; j++)
			{
				if (eph[i][j].valid)
				{
					//对于GPS来说是3600s
					//如果满足星历toe与传入的toe之间相差3600以内，选择该组星历计算卫星位置
					if (abs(eph[i][j].sow - toe_Input) <= 3600.0)
					{
						//时间混乱，需要搞清楚。
						//这组星历eph[i][j]是满足要求的，用它计算
						//卫星钟差
						MJD mjd_transmit;
						mjd_transmit.sod = obst.mjd.sod;
						mjd_transmit.Days = obst.mjd.Days;
						double deltaT = 0.0;						//t1 - t0
						double SatClk0 = eph[i][j].a0;
						double SatClk = eph[i][j].a0;
						double deltaTr = CaldeltaTr(mjd_transmit, eph[i][j]);
						do {
							SatClk0 = SatClk;
							mjd_transmit.sod = obst.mjd.sod - P_if / CLIGHT - SatClk0;
							deltaTr = CaldeltaTr(mjd_transmit, eph[i][j]);
							deltaT = mjd_transmit.ToGPST().sow - eph[i][j].sow; //t1 - toe
							SatClk = eph[i][j].a0 + eph[i][j].a1*deltaT + eph[i][j].a2 * deltaT*deltaT + deltaTr;

						} while (abs(SatClk - SatClk0) > 1e-9);
						sat[i].transmitTime = mjd_transmit;
						sat[i].prn = i;
						sat[i].ttr = obst.mjd.sod - mjd_transmit.sod;
						sat[i].satclk = SatClk;
						//计算信号发射时刻的卫星位置
						calSingleSatPV(mjd_transmit, eph[i][j], sat[i]);
						sat[i].valid = true;
					}
				}
			}
		}
	}
}


double ComputeE(double Mt, double e)
{
	double Et0 = Mt;
	double Et;
	double deltaEt;
	do {
		Et = Mt + e*sin(Et0);
		deltaEt = abs(Et - Et0);
		Et0 = Et;
	} while (!(deltaEt < 1e-15));
	return Et;
}

/***********************************************************
 * func:   北斗卫星位置计算，给定某颗卫星的一组星历，算出指定时间的卫星位置
 * @paras:  [in]  inputT,MJD格式，指定的时间
 *		    [in]  eph,给定卫星的一组星历
 *          [out] sat,保存计算出的卫星位置，置标志位为真		
 * @return: void
***********************************************************/
void calSingleSatPVBDS(MJD imputT, const eph_t& eph, sat_t& sat)
{
	// 卫星位置根据啥时间来算，toe？

	//需要将传进来的mjd转成sow 周内秒 转成BDT的周内秒
	double t = imputT.ToGPST().sow;
	double toe = eph.sow;
	//平均角速度n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(MU_BDS) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//计算观测瞬间 t 时刻卫星的平近点角
	//这里的t为周内秒，和参考时刻不同一周时需要使t-toe的绝对值在302400间
	double Mt = eph.M0 + n*(t - toe);
	//判断Mt是否在0-2pi之间
	while (Mt > 2 * PI_BDS)
	{
		Mt = Mt - 2 * PI_BDS;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI_BDS;
	}

	//计算偏近点角
	double Et = ComputeE(Mt, eph.e);

	//计算真近点角
	double upper = sin(Et)*sqrt(1 - eph.e*eph.e);
	double under = cos(Et) - eph.e;
	//double ft = atan(upper / under);
	double ft = atan2(upper, under);

	//计算升交距角（未经改正的）
	double ut_pie = eph.omega + ft;

	//计算卫星向径
	double rt_pie = eph.sqrtA*eph.sqrtA*(1 - eph.e*cos(Et));

	//计算摄动改正项
	double delta_ut = eph.cuc*cos(2 * ut_pie) + eph.cus*sin(2 * ut_pie);
	double delta_rt = eph.crc*cos(2 * ut_pie) + eph.crs*sin(2 * ut_pie);
	double delta_it = eph.cic*cos(2 * ut_pie) + eph.cis*sin(2 * ut_pie);

	//进行摄动改正
	double ut = ut_pie + delta_ut;
	double rt = rt_pie + delta_rt;
	double it = eph.i0 + eph.idot*(t - toe) + delta_it;

	//计算卫星在轨道平面坐标系中的位置
	double xt = rt*cos(ut);
	double yt = rt*sin(ut);

	//计算观测瞬间升交点经度
	//double Lt = eph.omega0 + (eph.omegaDot - OMGE_BDS)*t - eph.omegaDot*toe;
	double Lt = eph.omega0 + (eph.omegaDot - OMGE_BDS)*(t - toe) - OMGE_BDS*toe;
	//计算卫星在地固坐标系下的坐标
	sat.pos[0] = xt*cos(Lt) - yt*cos(it)*sin(Lt);
	sat.pos[1] = xt*sin(Lt) + yt*cos(it)*cos(Lt);
	sat.pos[2] = yt*sin(it);

	//计算钟差
	//sat.satclk = eph.a0 + eph.a1*(t - toe) + eph.a2*(t - toe)*(t - toe);
}

/***********************************************************
* func:   北斗卫星位置计算，给定某颗卫星的一组星历，算出指定时间的卫星位置
* @paras:  [in]  inputT,MJD格式，指定的时间
*		   [in]  eph,给定卫星的一组星历
*          [out] sat,保存计算出的卫星位置，置标志位为真
* @return: void
***********************************************************/
void calSingleSatPVBDSGEO(MJD imputT, const eph_t& eph, sat_t& sat)
{
	//需要将传进来的mjd转成sow 周内秒 转成BDT的周内秒
	double t = imputT.ToGPST().sow ;//如果传进来的是BDT的MJD，那么t也是BDT
	double toe = eph.sow;
	//平均角速度n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(MU_BDS) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//计算观测瞬间 t 时刻卫星的平近点角
	//这里的t为周内秒，和参考时刻不同一周时需要使t-toe的绝对值在302400间
	double Mt = eph.M0 + n*(t - toe);
	//判断Mt是否在0-2pi之间
	while (Mt > 2 * PI_BDS)
	{
		Mt = Mt - 2 * PI_BDS;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI_BDS;
	}

	//计算偏近点角
	double Et = ComputeE(Mt, eph.e);

	//计算真近点角
	double upper = sin(Et)*sqrt(1 - eph.e*eph.e);
	double under = cos(Et) - eph.e;
	//double ft = atan(upper / under);
	double ft = atan2(upper, under);

	//计算升交距角（未经改正的）
	double ut_pie = eph.omega + ft;

	//计算卫星向径
	double rt_pie = eph.sqrtA*eph.sqrtA*(1 - eph.e*cos(Et));

	//计算摄动改正项
	double delta_ut = eph.cuc*cos(2 * ut_pie) + eph.cus*sin(2 * ut_pie);
	double delta_rt = eph.crc*cos(2 * ut_pie) + eph.crs*sin(2 * ut_pie);
	double delta_it = eph.cic*cos(2 * ut_pie) + eph.cis*sin(2 * ut_pie);

	//进行摄动改正
	double ut = ut_pie + delta_ut;
	double rt = rt_pie + delta_rt;
	double it = eph.i0 + eph.idot*(t - toe) + delta_it;

	//计算卫星在轨道平面坐标系中的位置
	double xt = rt*cos(ut);
	double yt = rt*sin(ut);

	//计算观测瞬间升交点经度 (惯性系)
	//double Lt = eph.omega0 + eph.omegaDot*t - eph.omegaDot*toe;
	double Lt = eph.omega0 + eph.omegaDot*(t - toe) - OMGE_BDS*toe;

	//计算GEO卫星在自定义坐标系下的坐标
	double x_gk = xt*cos(Lt) - yt*cos(it)*sin(Lt);
	double y_gk = xt*sin(Lt) + yt*cos(it)*cos(Lt);
	double z_gk = yt*sin(it);

	//计算GEO卫星在BDCS坐标系中的坐标
	double phi_z = (t-toe)*OMGE_BDS;
	double phi_x = -5.0 * PI_BDS / 180.0;

	
	sat.pos[0] = cos(phi_z)*x_gk + sin(phi_z)*cos(phi_x)*y_gk + sin(phi_z)*sin(phi_x)*z_gk;
	sat.pos[1] = -sin(phi_z)*x_gk + cos(phi_z)*cos(phi_x)*y_gk + cos(phi_z)*sin(phi_x)*z_gk;
	sat.pos[2] = -sin(phi_x)*y_gk + cos(phi_x)*z_gk;

	//计算钟差
	//sat.satclk = eph.a0 + eph.a1*(t - toe) + eph.a2*(t - toe)*(t - toe);
}

/***********************************************************
* func:    北斗.计算卫星钟差相对论改正项deltaTr
* @paras:  [in]  inputT,MJD格式，指定的时间
*		   [in]  eph,给定卫星的一组星历
* @return: 卫星钟差相对论改正项deltaTr
***********************************************************/
double CaldeltaTrBDS(MJD imputT, const eph_t& eph)
{
	//需要将传进来的mjd转成sow 周内秒
	double t = imputT.ToGPST().sow;
	double toe = eph.sow;
	//平均角速度n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(MU_BDS) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//计算观测瞬间 t 时刻卫星的平近点角
	//这里的t为周内秒，和参考时刻不同一周时需要使t-toe的绝对值在302400间
	double Mt = eph.M0 + n*(t - toe);
	//判断Mt是否在0-2pi之间
	while (Mt > 2 * PI_BDS)
	{
		Mt = Mt - 2 * PI_BDS;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI_BDS;
	}

	//计算偏近点角
	double Et = ComputeE(Mt, eph.e);

	double deltaTr = -2 * sqrt(eph.sqrtA*eph.sqrtA*MU_BDS) * eph.e*sin(Et) / CLIGHT / CLIGHT;
	return deltaTr;
}

/***********************************************************
* func:    北斗.给定参考时间，批量计算卫星位置
* @paras:  [in]  inputT,MJD格式，指定的时间
*		   [in]  eph,卫星的星历
*		   [out] 存放计算出的卫星位置，同时标志位置为true
* @return: void
***********************************************************/
// input放的应该是obs.mjd，所以sat存放的是那个历元所有卫星的坐标
void CalBatchSatPVBDS(MJD inputT, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM])
{
	//观测值的toe，转成BDT的周内秒
	double toe_Input = inputT.ToGPST().sow;

	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		//寻找最近的星历，对某颗卫星而言
		for (int j = 0; j < EPOCHNUM; j++)
		{
			if (eph[i][j].valid)
			{
				//对于BDS来说是1800s
				//如果满足星历toe与传入的toe之间相差1800以内，选择该组星历计算卫星位置
				if (abs(eph[i][j].sow - toe_Input) <= 1800.0)
				{
					calSingleSatPV(inputT, eph[i][j], sat[i]);
					sat[i].valid = true;
				}
			}
		}
	}
}


/***********************************************************
* func:    北斗 迭代计算卫星信号发射时刻的卫星位置
* @paras:  [in]  obst,观测时刻及观测值，时间系统是GPST
*		   [in]  eph,卫星的星历
*		   [out] 存放计算出的卫星位置，同时标志位置为true
* @return: void
* Caution：观测文件里是GPST，计算北斗卫星位置的时候要转成BDT，即在
*		   GPST基础上加上14s
***********************************************************/
void calSatTransPVBDS(obs_t& obst, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM])
{
	//根据接收到的时间来确定发射时间

	//梅开二度，需要将所有标志为置成false
	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		sat[i].valid = false;
	}

	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		if (obst.obsd[i].valid)		//判断他有效
		{
			//这里可以用双频消电离层伪距组合 Page107 对于GPS和BDS可能不一样吧
			//double P_if = 2.54573* obst.obsd[i].P[0] - 1.54573* obst.obsd[i].P[1];
			double v1 = SQR(BDSB1FREQUENCE) / (SQR(BDSB1FREQUENCE) - SQR(BDSB3FREQUENCE));
			double v2 = SQR(BDSB3FREQUENCE) / (SQR(BDSB1FREQUENCE) - SQR(BDSB3FREQUENCE));
			double P_if = v1 * obst.obsd[i].P[0] - v2 * obst.obsd[i].P[1];
			double toe_Input = obst.mjd.ToGPST().sow - 14; //转成BDT的sow

			//寻找最近的星历，对某颗卫星而言
			for (int j = 0; j < EPOCHNUM; j++)
			{
				if (eph[i][j].valid)
				{
					//对于GPS来说是3600s
					//如果满足星历toe与传入的toe之间相差3600以内，选择该组星历计算卫星位置
					if (abs(eph[i][j].sow - toe_Input) <= 1800.0)
					{
						//时间混乱，需要搞清楚。
						//这组星历eph[i][j]是满足要求的，用它计算
						//卫星钟差
						MJD mjd_transmit;
						mjd_transmit.sod = obst.mjd.sod - 14; //转成BDT的sod
						mjd_transmit.Days = obst.mjd.Days;
						double deltaT = 0.0;						//t1 - t0
						double SatClk0 = eph[i][j].a0;
						double SatClk = eph[i][j].a0;
						double deltaTr = CaldeltaTrBDS(mjd_transmit, eph[i][j]);
						do {
							SatClk0 = SatClk;
							mjd_transmit.sod = obst.mjd.sod - 14 - P_if / CLIGHT - SatClk0;
							deltaTr = CaldeltaTrBDS(mjd_transmit, eph[i][j]);
							deltaT = mjd_transmit.ToGPST().sow  - eph[i][j].sow; //t1 - toe 转成BDT
							SatClk = eph[i][j].a0 + eph[i][j].a1*deltaT + eph[i][j].a2 * deltaT*deltaT + deltaTr;

						} while (abs(SatClk - SatClk0) > 1e-9);
						sat[i].transmitTime = mjd_transmit;//BDT
						sat[i].prn = i;
						sat[i].ttr = obst.mjd.sod - 14 - mjd_transmit.sod;
						sat[i].satclk = SatClk;
						//计算信号发射时刻的卫星位置 TO CHECK应该没问题
						CommonTime caaa = mjd_transmit.ToCommonTime();
						if (i+1 < 6 + 32 || i+1>58 + 32)
							calSingleSatPVBDSGEO(mjd_transmit, eph[i][j], sat[i]);
						else
							calSingleSatPVBDS(mjd_transmit, eph[i][j], sat[i]);
						sat[i].valid = true;
					}
				}
			}
		}
	}
}