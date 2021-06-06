#include "calSatPV.h"

/************�����ڲ�ʹ�õ�С����*******************/
double ComputeE(double Mt, double e);	//��������Et

//����Ĭ�ϴ���������ĳ�����ǵ�һ����Ч����
void calSingleSatPV(MJD imputT, const eph_t& eph, sat_t& sat)
{
	// ����λ�ø���ɶʱ�����㣬toe��

	//��Ҫ����������mjdת��sow ������
	double t = imputT.ToGPST().sow;
	double toe = eph.sow;
	//ƽ�����ٶ�n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(mu) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//����۲�˲�� t ʱ�����ǵ�ƽ�����
	//�����tΪ�����룬�Ͳο�ʱ�̲�ͬһ��ʱ��Ҫʹt-toe�ľ���ֵ��302400��
	double Mt = eph.M0 + n*(t - toe);
	//�ж�Mt�Ƿ���0-2pi֮��
	while (Mt > 2 * PI)
	{
		Mt = Mt - 2 * PI;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI;
	}

	//����ƫ�����
	double Et = ComputeE(Mt,eph.e);

	//����������
	double upper = sin(Et)*sqrt(1 -eph.e*eph.e);
	double under = cos(Et) - eph.e;
	//double ft = atan(upper / under);
	double ft = atan2(upper, under);

	//����������ǣ�δ�������ģ�
	double ut_pie = eph.omega + ft;

	//����������
	double rt_pie = eph.sqrtA*eph.sqrtA*(1 - eph.e*cos(Et));

	//�����㶯������
	double delta_ut = eph.cuc*cos(2 * ut_pie) + eph.cus*sin(2 * ut_pie);
	double delta_rt = eph.crc*cos(2 * ut_pie) + eph.crs*sin(2 * ut_pie);
	double delta_it = eph.cic*cos(2 * ut_pie) + eph.cis*sin(2 * ut_pie);

	//�����㶯����
	double ut = ut_pie + delta_ut;
	double rt = rt_pie + delta_rt;
	double it = eph.i0 + eph.idot*(t - toe) + delta_it;

	//���������ڹ��ƽ������ϵ�е�λ��
	double xt = rt*cos(ut);
	double yt = rt*sin(ut);

	//����۲�˲�������㾭��
	double Lt = eph.omega0 + (eph.omegaDot - OMGE)*t - eph.omegaDot*toe;

	//���������ڵع�����ϵ�µ�����
	sat.pos[0] = xt*cos(Lt) - yt*cos(it)*sin(Lt);
	sat.pos[1] = xt*sin(Lt) + yt*cos(it)*cos(Lt);
	sat.pos[2] = yt*sin(it);

	//�����Ӳ�
	//sat.satclk = eph.a0 + eph.a1*(t - toe) + eph.a2*(t - toe)*(t - toe);
	
}

// ����deltaTr
double CaldeltaTr(MJD imputT, const eph_t& eph)
{
	//��Ҫ����������mjdת��sow ������
	double t = imputT.ToGPST().sow;
	double toe = eph.sow;
	//ƽ�����ٶ�n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(mu) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//����۲�˲�� t ʱ�����ǵ�ƽ�����
	//�����tΪ�����룬�Ͳο�ʱ�̲�ͬһ��ʱ��Ҫʹt-toe�ľ���ֵ��302400��
	double Mt = eph.M0 + n*(t - toe);
	//�ж�Mt�Ƿ���0-2pi֮��
	while (Mt > 2 * PI)
	{
		Mt = Mt - 2 * PI;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI;
	}

	//����ƫ�����
	double Et = ComputeE(Mt, eph.e);

	double deltaTr = -2 * sqrt(eph.sqrtA*eph.sqrtA*mu) * eph.e*sin(Et) / CLIGHT / CLIGHT;
	return deltaTr;
}

// input�ŵ�Ӧ����obs.mjd������sat��ŵ����Ǹ���Ԫ�������ǵ�����
void CalBatchSatPV(MJD inputT,  eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM])
{	
	//�۲�ֵ��toe
	double toe_Input = inputT.ToGPST().sow;

	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		//Ѱ���������������ĳ�����Ƕ���
		for (int j = 0; j < EPOCHNUM; j++)
		{
			if (eph[i][j].valid)
			{
				//����GPS��˵��3600s
				//�����������toe�봫���toe֮�����3600���ڣ�ѡ�����������������λ��
				if (abs(eph[i][j].sow - toe_Input) <= 3600.0)
				{
					calSingleSatPV(inputT, eph[i][j], sat[i]);
					sat[i].valid = true;
				}
			}
		}
	}
}

//�������������źŷ���ʱ�̵�����λ��,�о�����Ҫα�����ֵ��
void calSatTransPV(obs_t& obst, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM])
{
	//���ݽ��յ���ʱ����ȷ������ʱ��

	//÷�����ȣ���Ҫ�����б�־Ϊ�ó�false
	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		sat[i].valid = false;
	}

	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		if (obst.obsd[i].valid)		//�ж�����Ч
		{
			
			//���������˫Ƶ�������α����� Page107 ����GPS��BDS���ܲ�һ����
			double P_if = 2.54573* obst.obsd[i].P[0] - 1.54573* obst.obsd[i].P[1];
			double toe_Input = obst.mjd.ToGPST().sow;

			//Ѱ���������������ĳ�����Ƕ���
			for (int j = 0; j < EPOCHNUM; j++)
			{
				if (eph[i][j].valid)
				{
					//����GPS��˵��3600s
					//�����������toe�봫���toe֮�����3600���ڣ�ѡ�����������������λ��
					if (abs(eph[i][j].sow - toe_Input) <= 3600.0)
					{
						//ʱ����ң���Ҫ�������
						//��������eph[i][j]������Ҫ��ģ���������
						//�����Ӳ�
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
						//�����źŷ���ʱ�̵�����λ��
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
 * func:   ��������λ�ü��㣬����ĳ�����ǵ�һ�����������ָ��ʱ�������λ��
 * @paras:  [in]  inputT,MJD��ʽ��ָ����ʱ��
 *		    [in]  eph,�������ǵ�һ������
 *          [out] sat,��������������λ�ã��ñ�־λΪ��		
 * @return: void
***********************************************************/
void calSingleSatPVBDS(MJD imputT, const eph_t& eph, sat_t& sat)
{
	// ����λ�ø���ɶʱ�����㣬toe��

	//��Ҫ����������mjdת��sow ������ ת��BDT��������
	double t = imputT.ToGPST().sow;
	double toe = eph.sow;
	//ƽ�����ٶ�n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(MU_BDS) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//����۲�˲�� t ʱ�����ǵ�ƽ�����
	//�����tΪ�����룬�Ͳο�ʱ�̲�ͬһ��ʱ��Ҫʹt-toe�ľ���ֵ��302400��
	double Mt = eph.M0 + n*(t - toe);
	//�ж�Mt�Ƿ���0-2pi֮��
	while (Mt > 2 * PI_BDS)
	{
		Mt = Mt - 2 * PI_BDS;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI_BDS;
	}

	//����ƫ�����
	double Et = ComputeE(Mt, eph.e);

	//����������
	double upper = sin(Et)*sqrt(1 - eph.e*eph.e);
	double under = cos(Et) - eph.e;
	//double ft = atan(upper / under);
	double ft = atan2(upper, under);

	//����������ǣ�δ�������ģ�
	double ut_pie = eph.omega + ft;

	//����������
	double rt_pie = eph.sqrtA*eph.sqrtA*(1 - eph.e*cos(Et));

	//�����㶯������
	double delta_ut = eph.cuc*cos(2 * ut_pie) + eph.cus*sin(2 * ut_pie);
	double delta_rt = eph.crc*cos(2 * ut_pie) + eph.crs*sin(2 * ut_pie);
	double delta_it = eph.cic*cos(2 * ut_pie) + eph.cis*sin(2 * ut_pie);

	//�����㶯����
	double ut = ut_pie + delta_ut;
	double rt = rt_pie + delta_rt;
	double it = eph.i0 + eph.idot*(t - toe) + delta_it;

	//���������ڹ��ƽ������ϵ�е�λ��
	double xt = rt*cos(ut);
	double yt = rt*sin(ut);

	//����۲�˲�������㾭��
	//double Lt = eph.omega0 + (eph.omegaDot - OMGE_BDS)*t - eph.omegaDot*toe;
	double Lt = eph.omega0 + (eph.omegaDot - OMGE_BDS)*(t - toe) - OMGE_BDS*toe;
	//���������ڵع�����ϵ�µ�����
	sat.pos[0] = xt*cos(Lt) - yt*cos(it)*sin(Lt);
	sat.pos[1] = xt*sin(Lt) + yt*cos(it)*cos(Lt);
	sat.pos[2] = yt*sin(it);

	//�����Ӳ�
	//sat.satclk = eph.a0 + eph.a1*(t - toe) + eph.a2*(t - toe)*(t - toe);
}

/***********************************************************
* func:   ��������λ�ü��㣬����ĳ�����ǵ�һ�����������ָ��ʱ�������λ��
* @paras:  [in]  inputT,MJD��ʽ��ָ����ʱ��
*		   [in]  eph,�������ǵ�һ������
*          [out] sat,��������������λ�ã��ñ�־λΪ��
* @return: void
***********************************************************/
void calSingleSatPVBDSGEO(MJD imputT, const eph_t& eph, sat_t& sat)
{
	//��Ҫ����������mjdת��sow ������ ת��BDT��������
	double t = imputT.ToGPST().sow ;//�������������BDT��MJD����ôtҲ��BDT
	double toe = eph.sow;
	//ƽ�����ٶ�n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(MU_BDS) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//����۲�˲�� t ʱ�����ǵ�ƽ�����
	//�����tΪ�����룬�Ͳο�ʱ�̲�ͬһ��ʱ��Ҫʹt-toe�ľ���ֵ��302400��
	double Mt = eph.M0 + n*(t - toe);
	//�ж�Mt�Ƿ���0-2pi֮��
	while (Mt > 2 * PI_BDS)
	{
		Mt = Mt - 2 * PI_BDS;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI_BDS;
	}

	//����ƫ�����
	double Et = ComputeE(Mt, eph.e);

	//����������
	double upper = sin(Et)*sqrt(1 - eph.e*eph.e);
	double under = cos(Et) - eph.e;
	//double ft = atan(upper / under);
	double ft = atan2(upper, under);

	//����������ǣ�δ�������ģ�
	double ut_pie = eph.omega + ft;

	//����������
	double rt_pie = eph.sqrtA*eph.sqrtA*(1 - eph.e*cos(Et));

	//�����㶯������
	double delta_ut = eph.cuc*cos(2 * ut_pie) + eph.cus*sin(2 * ut_pie);
	double delta_rt = eph.crc*cos(2 * ut_pie) + eph.crs*sin(2 * ut_pie);
	double delta_it = eph.cic*cos(2 * ut_pie) + eph.cis*sin(2 * ut_pie);

	//�����㶯����
	double ut = ut_pie + delta_ut;
	double rt = rt_pie + delta_rt;
	double it = eph.i0 + eph.idot*(t - toe) + delta_it;

	//���������ڹ��ƽ������ϵ�е�λ��
	double xt = rt*cos(ut);
	double yt = rt*sin(ut);

	//����۲�˲�������㾭�� (����ϵ)
	//double Lt = eph.omega0 + eph.omegaDot*t - eph.omegaDot*toe;
	double Lt = eph.omega0 + eph.omegaDot*(t - toe) - OMGE_BDS*toe;

	//����GEO�������Զ�������ϵ�µ�����
	double x_gk = xt*cos(Lt) - yt*cos(it)*sin(Lt);
	double y_gk = xt*sin(Lt) + yt*cos(it)*cos(Lt);
	double z_gk = yt*sin(it);

	//����GEO������BDCS����ϵ�е�����
	double phi_z = (t-toe)*OMGE_BDS;
	double phi_x = -5.0 * PI_BDS / 180.0;

	
	sat.pos[0] = cos(phi_z)*x_gk + sin(phi_z)*cos(phi_x)*y_gk + sin(phi_z)*sin(phi_x)*z_gk;
	sat.pos[1] = -sin(phi_z)*x_gk + cos(phi_z)*cos(phi_x)*y_gk + cos(phi_z)*sin(phi_x)*z_gk;
	sat.pos[2] = -sin(phi_x)*y_gk + cos(phi_x)*z_gk;

	//�����Ӳ�
	//sat.satclk = eph.a0 + eph.a1*(t - toe) + eph.a2*(t - toe)*(t - toe);
}

/***********************************************************
* func:    ����.���������Ӳ�����۸�����deltaTr
* @paras:  [in]  inputT,MJD��ʽ��ָ����ʱ��
*		   [in]  eph,�������ǵ�һ������
* @return: �����Ӳ�����۸�����deltaTr
***********************************************************/
double CaldeltaTrBDS(MJD imputT, const eph_t& eph)
{
	//��Ҫ����������mjdת��sow ������
	double t = imputT.ToGPST().sow;
	double toe = eph.sow;
	//ƽ�����ٶ�n0
	double sqrt_A3 = eph.sqrtA * eph.sqrtA * eph.sqrtA;
	double n0 = sqrt(MU_BDS) / (sqrt_A3);
	double n = n0 + eph.delta_n;

	//����۲�˲�� t ʱ�����ǵ�ƽ�����
	//�����tΪ�����룬�Ͳο�ʱ�̲�ͬһ��ʱ��Ҫʹt-toe�ľ���ֵ��302400��
	double Mt = eph.M0 + n*(t - toe);
	//�ж�Mt�Ƿ���0-2pi֮��
	while (Mt > 2 * PI_BDS)
	{
		Mt = Mt - 2 * PI_BDS;
	}
	while (Mt < 0)
	{
		Mt = Mt + 2 * PI_BDS;
	}

	//����ƫ�����
	double Et = ComputeE(Mt, eph.e);

	double deltaTr = -2 * sqrt(eph.sqrtA*eph.sqrtA*MU_BDS) * eph.e*sin(Et) / CLIGHT / CLIGHT;
	return deltaTr;
}

/***********************************************************
* func:    ����.�����ο�ʱ�䣬������������λ��
* @paras:  [in]  inputT,MJD��ʽ��ָ����ʱ��
*		   [in]  eph,���ǵ�����
*		   [out] ��ż����������λ�ã�ͬʱ��־λ��Ϊtrue
* @return: void
***********************************************************/
// input�ŵ�Ӧ����obs.mjd������sat��ŵ����Ǹ���Ԫ�������ǵ�����
void CalBatchSatPVBDS(MJD inputT, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM])
{
	//�۲�ֵ��toe��ת��BDT��������
	double toe_Input = inputT.ToGPST().sow;

	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		//Ѱ���������������ĳ�����Ƕ���
		for (int j = 0; j < EPOCHNUM; j++)
		{
			if (eph[i][j].valid)
			{
				//����BDS��˵��1800s
				//�����������toe�봫���toe֮�����1800���ڣ�ѡ�����������������λ��
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
* func:    ���� �������������źŷ���ʱ�̵�����λ��
* @paras:  [in]  obst,�۲�ʱ�̼��۲�ֵ��ʱ��ϵͳ��GPST
*		   [in]  eph,���ǵ�����
*		   [out] ��ż����������λ�ã�ͬʱ��־λ��Ϊtrue
* @return: void
* Caution���۲��ļ�����GPST�����㱱������λ�õ�ʱ��Ҫת��BDT������
*		   GPST�����ϼ���14s
***********************************************************/
void calSatTransPVBDS(obs_t& obst, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM])
{
	//���ݽ��յ���ʱ����ȷ������ʱ��

	//÷�����ȣ���Ҫ�����б�־Ϊ�ó�false
	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		sat[i].valid = false;
	}

	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		if (obst.obsd[i].valid)		//�ж�����Ч
		{
			//���������˫Ƶ�������α����� Page107 ����GPS��BDS���ܲ�һ����
			//double P_if = 2.54573* obst.obsd[i].P[0] - 1.54573* obst.obsd[i].P[1];
			double v1 = SQR(BDSB1FREQUENCE) / (SQR(BDSB1FREQUENCE) - SQR(BDSB3FREQUENCE));
			double v2 = SQR(BDSB3FREQUENCE) / (SQR(BDSB1FREQUENCE) - SQR(BDSB3FREQUENCE));
			double P_if = v1 * obst.obsd[i].P[0] - v2 * obst.obsd[i].P[1];
			double toe_Input = obst.mjd.ToGPST().sow - 14; //ת��BDT��sow

			//Ѱ���������������ĳ�����Ƕ���
			for (int j = 0; j < EPOCHNUM; j++)
			{
				if (eph[i][j].valid)
				{
					//����GPS��˵��3600s
					//�����������toe�봫���toe֮�����3600���ڣ�ѡ�����������������λ��
					if (abs(eph[i][j].sow - toe_Input) <= 1800.0)
					{
						//ʱ����ң���Ҫ�������
						//��������eph[i][j]������Ҫ��ģ���������
						//�����Ӳ�
						MJD mjd_transmit;
						mjd_transmit.sod = obst.mjd.sod - 14; //ת��BDT��sod
						mjd_transmit.Days = obst.mjd.Days;
						double deltaT = 0.0;						//t1 - t0
						double SatClk0 = eph[i][j].a0;
						double SatClk = eph[i][j].a0;
						double deltaTr = CaldeltaTrBDS(mjd_transmit, eph[i][j]);
						do {
							SatClk0 = SatClk;
							mjd_transmit.sod = obst.mjd.sod - 14 - P_if / CLIGHT - SatClk0;
							deltaTr = CaldeltaTrBDS(mjd_transmit, eph[i][j]);
							deltaT = mjd_transmit.ToGPST().sow  - eph[i][j].sow; //t1 - toe ת��BDT
							SatClk = eph[i][j].a0 + eph[i][j].a1*deltaT + eph[i][j].a2 * deltaT*deltaT + deltaTr;

						} while (abs(SatClk - SatClk0) > 1e-9);
						sat[i].transmitTime = mjd_transmit;//BDT
						sat[i].prn = i;
						sat[i].ttr = obst.mjd.sod - 14 - mjd_transmit.sod;
						sat[i].satclk = SatClk;
						//�����źŷ���ʱ�̵�����λ�� TO CHECKӦ��û����
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