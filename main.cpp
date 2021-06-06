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
	//结果文件输出
	FILE *pfile = fopen(kFloatResOutFilepath, "w");
	FILE *pfilefix = fopen(kFixResOutFilepath, "w");

	/**************结构体定义,基站和流动站观测值、头、星历，卫星位置sat，计算结果sol*************/
	obs_t obs_base, obs_rover;
	obs_header obsHeader_base, obsHeader_rover;
	nav_header navHeader;
	//下面使用new的方式减少堆栈的内存，避免堆栈溢出
	eph_TArray *eph = new eph_t[MAXSATNUM][EPOCHNUM];
	sat_t *sat_base = new sat_t[MAXSATNUM];
	sat_t *sat_rover = new sat_t[MAXSATNUM];

	sol_t sol_base;						//SPP每个历元的解算结果
	sol_t sol_rover;					//流动站解算结果
	/*===================================================================================*/

	
	/******************************读取导航电文Nav存入eph中*************************************/
	//novatel 0baseline filepath: NovAtel719Data1_20210406\\brdm0960.21p
	fstream fin(kEphFilepath, ios::in);
	readNavHeader(fin, navHeader);
	readNav(fin, eph);
	fin.close();
	/**===================================================================================***/

	fstream fs_obs_base(kBaseFilepath, ios::in);	//只读方式打开
	fstream fs_obs_rover(kRoverFilepath, ios::in);	//只读方式打开
	if (!fs_obs_base.is_open() || !fs_obs_rover.is_open())						//判断文件正常打开,To Improve
	{	
		printf("ERROR");					// to do 返回错误信息
	}

	/***************读取观测值头文件***********************************************/
	readObsHeader(fs_obs_base, obsHeader_base);
	readObsHeader(fs_obs_rover, obsHeader_rover);

	//使用O文件里的概略坐标，但不是所有文件都存有，因此暂未使用
	//sol_base.pos[0] = obsHeader_base.approxpos_rcv_xyz[0];
	//sol_base.pos[1] = obsHeader_base.approxpos_rcv_xyz[1];
	//sol_base.pos[2] = obsHeader_base.approxpos_rcv_xyz[2];
	//sol_rover.pos[0] = obsHeader_rover.approxpos_rcv_xyz[0];
	//sol_rover.pos[1] = obsHeader_rover.approxpos_rcv_xyz[1];
	//sol_rover.pos[2] = obsHeader_rover.approxpos_rcv_xyz[2];

	//精确坐标
	sol_base.pos[0] = kPreciseBaseCoor[0];
	sol_base.pos[1] = kPreciseBaseCoor[1];
	sol_base.pos[2] = kPreciseBaseCoor[2];
	sol_rover.pos[0] = kPreciseRoverCoor[0];
	sol_rover.pos[1] = kPreciseRoverCoor[1];
	sol_rover.pos[2] = kPreciseRoverCoor[2];

	/****************************************************************************/

	/*浮点解和固定解的输出*/
	RtkFloat rtk_float;
	RtkFix rtk_fix;

	//可以输出，但是有点问题，因为判断条件不太对
	//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@貌似结束不了循环？？？？ TO FIXED

	double timeGap = 0;
	const double requireGap = 0.01;
	bool jumpLoop = false;	//跳出循环的条件
	while (!fs_obs_base.eof() && !fs_obs_rover.eof()) //时间同步，要求是0.01s
	{
		/*========================时间同步=============================================*/
		readObs(fs_obs_base, obs_base);
		readObs(fs_obs_rover, obs_rover);
		timeGap = (int)(obs_base.ct.Minute * 60 + obs_base.ct.Second - (obs_rover.ct.Minute * 60 + obs_rover.ct.Second));
		if (timeGap < 0) //说明基站在前，流动站在后，基站慢了
		{
			//如果两者时间差大于0.01s，则继续读取直至同步。
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
		/*========================时间同步完成=============================================*/

		//计算卫星位置 应该没啥问题了
		calSatTransPV(obs_base, eph, sat_base);
		calSatTransPV(obs_rover, eph, sat_rover);
		calSatTransPVBDS(obs_base, eph, sat_base);
		calSatTransPVBDS(obs_rover, eph, sat_rover);
		//RTK函数模型构建
		
		double delta_z = 0.0;
		double sol_rover_pos2 = sol_rover.pos[2];
		//需要迭代计算浮点解，直至两次sol_base.pos差值小于阈值
		do
		{
			sol_rover_pos2 = sol_rover.pos[2];
			CreateRTKModel(obs_base, obs_rover, sat_base, sat_rover, sol_base, sol_rover,rtk_float);
			delta_z = abs(sol_rover.pos[2] - sol_rover_pos2);
		} while (delta_z > 1e-4);
		fprintf(pfile, " % 7.4f % 7.4f % 7.4f\n", sol_rover.pos[0]-sol_base.pos[0], sol_rover.pos[1] - sol_base.pos[1], sol_rover.pos[2] - sol_base.pos[2]);
		FixAmbiguity(rtk_float, rtk_fix);
		//输出固定解
		fprintf(pfilefix, " % 7.4f % 7.4f % 7.4f %5f\n", rtk_fix._baseline_fix[0], rtk_fix._baseline_fix[1], rtk_fix._baseline_fix[2], rtk_fix.ratio);
		//printf(" % 7.4f % 7.4f % 7.4f %5f\n", rtk_fix._baseline_fix[0], rtk_fix._baseline_fix[1], rtk_fix._baseline_fix[2], rtk_fix.ratio);

		//迭代计算测站位置 SPP  To do
		
/*		double satX = sol_base.pos[0];
		do {
			LeastSquare(obs_base, sat_base, sol_base);
		} while (abs(satX - sol_base.pos[0]) > 1.0);
		
		printf("%f  %f  %f\n", sol_base.pos[0], sol_base.pos[1], sol_base.pos[2]);	*/	
	}

	/***********************************关闭文件*************************************/
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