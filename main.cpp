#include <iostream>

#include "calSatPV.h"
#include "config.h"
#include "global.h"
#include "readRinex.h"
#include "rinex.h"
#include "RTKSolution.h"
#include "SPPSolution.h"
#include "TimeCoor.h"

using namespace std;

eph_t a[EPOCHNUM];
typedef  decltype(a) eph_TArray;

int main()
{
	//����ļ����
	FILE *pfile = fopen(kFloatResOutFilepath, "w");
	FILE *pfilefix = fopen(kFixResOutFilepath, "w");

	/**************�ṹ�嶨��,��վ������վ�۲�ֵ��ͷ������������λ��sat��������sol*************/
	obs_t obs_base, obs_rover;
	obs_header obsHeader_base, obsHeader_rover;
	nav_header navHeader;
	//����ʹ��new�ķ�ʽ���ٶ�ջ���ڴ棬�����ջ���
	eph_TArray *eph = new eph_t[MAXSATNUM][EPOCHNUM];
	sat_t *sat_base = new sat_t[MAXSATNUM];
	sat_t *sat_rover = new sat_t[MAXSATNUM];

	sol_t sol_base;						//SPPÿ����Ԫ�Ľ�����
	sol_t sol_rover;					//����վ������
	/*===================================================================================*/

	
	/******************************��ȡ��������Nav����eph��*************************************/
	//novatel 0baseline filepath: NovAtel719Data1_20210406\\brdm0960.21p
	fstream fin(kEphFilepath, ios::in);
	readNavHeader(fin, navHeader);
	readNav(fin, eph);
	fin.close();
	/**===================================================================================***/

	fstream fs_obs_base(kBaseFilepath, ios::in);	//ֻ����ʽ��
	fstream fs_obs_rover(kRoverFilepath, ios::in);	//ֻ����ʽ��
	if (!fs_obs_base.is_open() || !fs_obs_rover.is_open())						//�ж��ļ�������,To Improve
	{	
		printf("ERROR");					// to do ���ش�����Ϣ
	}

	/***************��ȡ�۲�ֵͷ�ļ�***********************************************/
	readObsHeader(fs_obs_base, obsHeader_base);
	readObsHeader(fs_obs_rover, obsHeader_rover);

	//ʹ��O�ļ���ĸ������꣬�����������ļ������У������δʹ��
	//sol_base.pos[0] = obsHeader_base.approxpos_rcv_xyz[0];
	//sol_base.pos[1] = obsHeader_base.approxpos_rcv_xyz[1];
	//sol_base.pos[2] = obsHeader_base.approxpos_rcv_xyz[2];
	//sol_rover.pos[0] = obsHeader_rover.approxpos_rcv_xyz[0];
	//sol_rover.pos[1] = obsHeader_rover.approxpos_rcv_xyz[1];
	//sol_rover.pos[2] = obsHeader_rover.approxpos_rcv_xyz[2];

	//��ȷ����
	sol_base.pos[0] = kPreciseBaseCoor[0];
	sol_base.pos[1] = kPreciseBaseCoor[1];
	sol_base.pos[2] = kPreciseBaseCoor[2];
	sol_rover.pos[0] = kPreciseRoverCoor[0];
	sol_rover.pos[1] = kPreciseRoverCoor[1];
	sol_rover.pos[2] = kPreciseRoverCoor[2];

	/****************************************************************************/

	/*�����͹̶�������*/
	RtkFloat rtk_float;
	RtkFix rtk_fix;

	//��������������е����⣬��Ϊ�ж�������̫��
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ò�ƽ�������ѭ���������� TO FIXED

	double timeGap = 0;
	const double requireGap = 0.01;
	bool jumpLoop = false;	//����ѭ��������
	while (!fs_obs_base.eof() && !fs_obs_rover.eof()) //ʱ��ͬ����Ҫ����0.01s
	{
		/*========================ʱ��ͬ��=============================================*/
		readObs(fs_obs_base, obs_base);
		readObs(fs_obs_rover, obs_rover);
		timeGap = (int)(obs_base.ct.Minute * 60 + obs_base.ct.Second - (obs_rover.ct.Minute * 60 + obs_rover.ct.Second));
		if (timeGap < 0) //˵����վ��ǰ������վ�ں󣬻�վ����
		{
			//�������ʱ������0.01s���������ȡֱ��ͬ����
			for (; abs(timeGap) !=0 ;timeGap+=obsHeader_base.interval)
				if (!fs_obs_base.eof())	readObs(fs_obs_base, obs_base);
				else {
					jumpLoop = true;
					break;
				}
		}
		else if (timeGap > 0)
		{
			for (; abs(timeGap) != 0; timeGap-= obsHeader_base.interval)
				if (!fs_obs_rover.eof()) readObs(fs_obs_rover, obs_rover);
				else {
					jumpLoop = true;
					break;
				}
		}
		else;
		if (jumpLoop) break;
		/*========================ʱ��ͬ�����=============================================*/

		//��������λ�� Ӧ��ûɶ������
		calSatTransPV(obs_base, eph, sat_base);
		calSatTransPV(obs_rover, eph, sat_rover);
		calSatTransPVBDS(obs_base, eph, sat_base);
		calSatTransPVBDS(obs_rover, eph, sat_rover);
		//RTK����ģ�͹���
		
		double delta_z = 0.0;
		double sol_rover_pos2 = sol_rover.pos[2];
		//��Ҫ�������㸡��⣬ֱ������sol_base.pos��ֵС����ֵ
		do
		{
			sol_rover_pos2 = sol_rover.pos[2];
			CreateRTKModel(obs_base, obs_rover, sat_base, sat_rover, sol_base, sol_rover,rtk_float);
			delta_z = abs(sol_rover.pos[2] - sol_rover_pos2);
		} while (delta_z > 1e-4);
		fprintf(pfile, " % 7.4f % 7.4f % 7.4f\n", sol_rover.pos[0]-sol_base.pos[0], sol_rover.pos[1] - sol_base.pos[1], sol_rover.pos[2] - sol_base.pos[2]);
		FixAmbiguity(rtk_float, rtk_fix);
		//����̶���
		fprintf(pfilefix, " % 7.4f % 7.4f % 7.4f %5f\n", rtk_fix._baseline_fix[0], rtk_fix._baseline_fix[1], rtk_fix._baseline_fix[2], rtk_fix.ratio);
		//printf(" % 7.4f % 7.4f % 7.4f %5f\n", rtk_fix._baseline_fix[0], rtk_fix._baseline_fix[1], rtk_fix._baseline_fix[2], rtk_fix.ratio);

		//���������վλ�� SPP  To do
		
/*		double satX = sol_base.pos[0];
		do {
			LeastSquare(obs_base, sat_base, sol_base);
		} while (abs(satX - sol_base.pos[0]) > 1.0);
		
		printf("%f  %f  %f\n", sol_base.pos[0], sol_base.pos[1], sol_base.pos[2]);	*/	
	}

	/***********************************�ر��ļ�*************************************/
	fclose(pfile);
	fclose(pfilefix);
	fs_obs_base.close();
	fs_obs_rover.close();
	delete[] eph;
	delete[] sat_base;
	delete[] sat_rover;
	/*******************************************************************************/

	return 0;
}