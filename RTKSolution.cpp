#include "RTKSolution.h"

//高度角法算sigma^2, E是高度角，单位rad
double CalSigma2(double E);

//返回B、L
void CreateRTKModel(obs_t& obs_base, obs_t& obs_rover, sat_t sat_base[MAXSATNUM],sat_t sat_rover[MAXSATNUM],sol_t& sol_base, sol_t& sol_rover,
					RtkFloat& rtk_float)
{
	/*GPS*************************************************************************************/
	//共视卫星选取
	int refSatIndex_GPS = 0, refSatIndex_BDS = 0;
	double refElev_GPS = 0.0, refElev_BDS = 0.0;

	int obsNum = 0;	
	int obsNum_GPS = 0, obsNum_BDS = 0;

	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		//如果两者不能同时满足，则跳过该卫星
		if (!(sat_rover[i].valid && sat_base[i].valid))
			continue;
		//计算卫星高度角
		obs_base.obsd[i].elev = RAD2DEG * CalElevation(sat_base[i].pos, sol_base.pos);
		obs_rover.obsd[i].elev = RAD2DEG *CalElevation(sat_rover[i].pos, sol_rover.pos);
		//如果卫星高度角小于15度，则丢弃
		if (obs_base.obsd[i].elev < kElevationMaskAngle)
		{
			obs_base.obsd[i].valid = false;
			obs_rover.obsd[i].valid = false;
			sat_rover[i].valid = false;
			sat_base[i].valid = false;
			continue;
		}
		//基准星选取
		if (obs_base.obsd[i].elev > refElev_GPS)
		{
			refElev_GPS = obs_base.obsd[i].elev;
			refSatIndex_GPS = i;
		}
		obsNum_GPS++;
	}
	obsNum = obsNum_GPS;
	//现在obsNum是GPS卫星可用数目,那么对于GPS来说共有obsnum-1个双差观测值

	//求北斗可用卫星数目
	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		//如果两者不能同时满足，则跳过该卫星
		if (!(sat_rover[i].valid && sat_base[i].valid))
			continue;
		//计算卫星高度角
		obs_base.obsd[i].elev = RAD2DEG * CalElevation(sat_base[i].pos, sol_base.pos);
		obs_rover.obsd[i].elev = RAD2DEG *CalElevation(sat_rover[i].pos, sol_rover.pos);
		//如果卫星高度角小于15度，则丢弃
		if (obs_base.obsd[i].elev < kElevationMaskAngle)
		{
			obs_base.obsd[i].valid = false;
			obs_rover.obsd[i].valid = false;
			sat_rover[i].valid = false;
			sat_base[i].valid = false;
			continue;
		}
		obsNum_BDS++;
		//北斗基准星选取
		//如果是GEO卫星，不参与选取工作
		if (i - 32 + 1 < 6 || i - 32 + 1 > 58)
			continue;
		if (obs_base.obsd[i].elev > refElev_BDS)
		{
			refElev_BDS = obs_base.obsd[i].elev;
			refSatIndex_BDS = i;
		}
	}
	obsNum = obsNum_BDS + obsNum_GPS;


	//用动态矩阵，这样构造的时候默认是全0阵
	//Matrix B(4 * obsNum_GPS - 4, 3 + 2 * obsNum_GPS - 2);
	//Matrix l(4 * obsNum_GPS - 4, 1);
	//Matrix P(4 * obsNum_GPS - 4, 4 * obsNum_GPS - 4);
	//Matrix D(4 * obsNum_GPS - 4, 4 * obsNum_GPS - 4);
	Matrix B(4 * obsNum - 8, 3 + 2 * obsNum - 4);
	Matrix l(4 * obsNum - 8, 1);
	Matrix P(4 * obsNum - 8, 4 * obsNum - 8);
	Matrix D(4 * obsNum - 8, 4 * obsNum - 8);


	//构建单差双差观测值,注意载波是以m为单位
	double diffP0ObsStation_ref = obs_rover.obsd[refSatIndex_GPS].P[0] - obs_base.obsd[refSatIndex_GPS].P[0];
	double diffP1ObsStation_ref = obs_rover.obsd[refSatIndex_GPS].P[1] - obs_base.obsd[refSatIndex_GPS].P[1];
	double diffL0ObsStation_ref = (obs_rover.obsd[refSatIndex_GPS].L[0] - obs_base.obsd[refSatIndex_GPS].L[0])*GPSL1LAMBDA;
	double diffL1ObsStation_ref = (obs_rover.obsd[refSatIndex_GPS].L[1] - obs_base.obsd[refSatIndex_GPS].L[1])*GPSL2LAMBDA;

	double diff_p0_bds_ref = obs_rover.obsd[refSatIndex_BDS].P[0] - obs_base.obsd[refSatIndex_BDS].P[0];
	double diff_p1_bds_ref = obs_rover.obsd[refSatIndex_BDS].P[1] - obs_base.obsd[refSatIndex_BDS].P[1];
	double diff_l0_bds_ref = (obs_rover.obsd[refSatIndex_BDS].L[0] - obs_base.obsd[refSatIndex_BDS].L[0])*BDSB1LAMBDA;
	double diff_l1_bds_ref = (obs_rover.obsd[refSatIndex_BDS].L[1] - obs_base.obsd[refSatIndex_BDS].L[1])*BDSB3LAMBDA;

	int iter = 0;
	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		if (sat_rover[i].valid && sat_base[i].valid && i != refSatIndex_GPS)
		{
			//站间单差
			double diffP0ObsStation = obs_rover.obsd[i].P[0] - obs_base.obsd[i].P[0];
			double diffP1ObsStation = obs_rover.obsd[i].P[1] - obs_base.obsd[i].P[1];
			double diffL0ObsStation = (obs_rover.obsd[i].L[0] - obs_base.obsd[i].L[0])*GPSL1LAMBDA;
			double diffL1ObsStation = (obs_rover.obsd[i].L[1] - obs_base.obsd[i].L[1])*GPSL2LAMBDA;
			//星间求差双差
			double doubleDiffP0 = diffP0ObsStation - diffP0ObsStation_ref;
			double doubleDiffP1 = diffP1ObsStation - diffP1ObsStation_ref;
			double doubleDiffL0 = diffL0ObsStation - diffL0ObsStation_ref;
			double doubleDiffL1 = diffL1ObsStation - diffL1ObsStation_ref;

			//双差的计算值c，为了o-c
			double distance_kB = CalDistance(sat_rover[i].pos, sol_rover.pos);
			double distance_jB = CalDistance(sat_rover[refSatIndex_GPS].pos, sol_rover.pos);
			double distance_kA = CalDistance(sat_base[i].pos, sol_base.pos);
			double distance_jA = CalDistance(sat_base[refSatIndex_GPS].pos, sol_base.pos);
			double doubleDiff_c = distance_kB - distance_jB - distance_kA + distance_jA;

			//计算B矩阵的系数 A是sol_base B是rover	j是参考星	k是流动星

			double rhoBk = CalDistance(sol_rover.pos, sat_rover[i].pos);
			double rhoBj = CalDistance(sol_rover.pos, sat_rover[refSatIndex_GPS].pos);
			//BDS和GPS融合算，只是参考星选取的不一样。
			double coeffXB = -(sat_rover[i].pos[0] - sol_rover.pos[0]) / rhoBk + (sat_rover[refSatIndex_GPS].pos[0] - sol_rover.pos[0]) / rhoBj;
			double coeffYB = -(sat_rover[i].pos[1] - sol_rover.pos[1]) / rhoBk + (sat_rover[refSatIndex_GPS].pos[1] - sol_rover.pos[1]) / rhoBj;
			double coeffZB = -(sat_rover[i].pos[2] - sol_rover.pos[2]) / rhoBk + (sat_rover[refSatIndex_GPS].pos[2] - sol_rover.pos[2]) / rhoBj;

			//先伪距再相位
			B(iter, 0) = coeffXB;
			B(iter, 1) = coeffYB;
			B(iter, 2) = coeffZB;
			B(iter + obsNum_GPS - 1, 0) = coeffXB;
			B(iter + obsNum_GPS - 1, 1) = coeffYB;
			B(iter + obsNum_GPS - 1, 2) = coeffZB;

			B(2* (obsNum_GPS-1) + iter, 0) = coeffXB;
			B(2 * (obsNum_GPS - 1) + iter, 1) = coeffYB;
			B(2 * (obsNum_GPS - 1) + iter, 2) = coeffZB;
			B(2 * (obsNum_GPS - 1) + iter, 3 + iter) = GPSL1LAMBDA;
			B(3 * (obsNum_GPS - 1) + iter, 0) = coeffXB;
			B(3 * (obsNum_GPS - 1) + iter, 1) = coeffYB;
			B(3 * (obsNum_GPS - 1) + iter, 2) = coeffZB;
			B(3 * (obsNum_GPS - 1) + iter, 3 + iter + obsNum_GPS - 1) = GPSL2LAMBDA;

			//用单频
			l(iter, 0) = doubleDiffP0 - doubleDiff_c;
			l(iter + obsNum_GPS - 1,0) = doubleDiffP1 - doubleDiff_c;
			l(2 * (obsNum_GPS-1) + iter, 0) = doubleDiffL0 - doubleDiff_c;
			l(3 * (obsNum_GPS - 1) + iter, 0) = doubleDiffL1 - doubleDiff_c;
			iter++;
			
		}
		//下面打算做北斗的。
	}


	//寻找一个偏移量,相当于把BDS的B矩阵整体偏移
	//BDS的
	iter = 0;
	int row_offset = 4 * obsNum_GPS - 4;
	int col_offset = 2 * obsNum_GPS - 2;
	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		if (sat_rover[i].valid && sat_base[i].valid && i != refSatIndex_BDS)
		{
			//站间单差
			double diff_p0_bds = obs_rover.obsd[i].P[0] - obs_base.obsd[i].P[0];
			double diff_p1_bds = obs_rover.obsd[i].P[1] - obs_base.obsd[i].P[1];
			double diff_l0_bds = (obs_rover.obsd[i].L[0] - obs_base.obsd[i].L[0])*BDSB1LAMBDA;
			double diff_l1_bds = (obs_rover.obsd[i].L[1] - obs_base.obsd[i].L[1])*BDSB3LAMBDA;
			//星间求差双差
			double double_diff_p0 = diff_p0_bds - diff_p0_bds_ref;
			double double_diff_p1 = diff_p1_bds - diff_p1_bds_ref;
			double double_diff_l0 = diff_l0_bds - diff_l0_bds_ref;
			double double_diff_l1 = diff_l1_bds - diff_l1_bds_ref;

			//双差的计算值c，为了o-c
			double distance_kB = CalDistance(sat_rover[i].pos, sol_rover.pos);
			double distance_jB = CalDistance(sat_rover[refSatIndex_BDS].pos, sol_rover.pos);
			double distance_kA = CalDistance(sat_base[i].pos, sol_base.pos);
			double distance_jA = CalDistance(sat_base[refSatIndex_BDS].pos, sol_base.pos);
			double double_diff_c = distance_kB - distance_jB - distance_kA + distance_jA;

			//计算B矩阵的系数 A是sol_base B是rover	j是参考星	k是流动星

			double rhoBk = CalDistance(sol_rover.pos, sat_rover[i].pos);
			double rhoBj = CalDistance(sol_rover.pos, sat_rover[refSatIndex_BDS].pos);
			//BDS和GPS融合算，只是参考星选取的不一样。
			double coeffXB = -(sat_rover[i].pos[0] - sol_rover.pos[0]) / rhoBk + (sat_rover[refSatIndex_BDS].pos[0] - sol_rover.pos[0]) / rhoBj;
			double coeffYB = -(sat_rover[i].pos[1] - sol_rover.pos[1]) / rhoBk + (sat_rover[refSatIndex_BDS].pos[1] - sol_rover.pos[1]) / rhoBj;
			double coeffZB = -(sat_rover[i].pos[2] - sol_rover.pos[2]) / rhoBk + (sat_rover[refSatIndex_BDS].pos[2] - sol_rover.pos[2]) / rhoBj;

			//先伪距再相位
			B(iter + row_offset, 0) = coeffXB;
			B(iter + row_offset, 1) = coeffYB;
			B(iter + row_offset, 2) = coeffZB;
			B(iter + obsNum_BDS - 1 + row_offset, 0) = coeffXB;
			B(iter + obsNum_BDS - 1 + row_offset, 1) = coeffYB;
			B(iter + obsNum_BDS - 1 + row_offset, 2) = coeffZB;

			B(2 * (obsNum_BDS - 1) + iter + row_offset, 0) = coeffXB;
			B(2 * (obsNum_BDS - 1) + iter + row_offset, 1) = coeffYB;
			B(2 * (obsNum_BDS - 1) + iter + row_offset, 2) = coeffZB;
			B(2 * (obsNum_BDS - 1) + iter + row_offset, 3 + iter + col_offset) = BDSB1LAMBDA;
			B(3 * (obsNum_BDS - 1) + iter + row_offset, 0) = coeffXB;
			B(3 * (obsNum_BDS - 1) + iter + row_offset, 1) = coeffYB;
			B(3 * (obsNum_BDS - 1) + iter + row_offset, 2) = coeffZB;
			B(3 * (obsNum_BDS - 1) + iter + row_offset, 3 + iter + obsNum_BDS - 1 + col_offset) = BDSB3LAMBDA;

			//用单频
			l(iter + row_offset, 0) = double_diff_p0 - double_diff_c;
			l(iter + obsNum_BDS - 1 + row_offset, 0) = double_diff_p1 - double_diff_c;
			l(2 * (obsNum_BDS - 1) + iter + row_offset, 0) = double_diff_l0 - double_diff_c;
			l(3 * (obsNum_BDS - 1) + iter + row_offset, 0) = double_diff_l1 - double_diff_c;
			iter++;
		}
		//下面打算做北斗的。
	}

	//随机模型

	//参考星单差观测值的方差
	double multiplePL = 10000.0;	//伪距方差相对于相位的倍数
	double Q_S_ref_L = CalSigma2(obs_base.obsd[refSatIndex_GPS].elev*DEG2RAD) +
		CalSigma2(obs_rover.obsd[refSatIndex_GPS].elev*DEG2RAD);
	double Q_S_ref_P = Q_S_ref_L * multiplePL;

	int doubleDiffNum_GPS = obsNum_GPS - 1;
	int row = 0;
	for (int i = 0; i < MAXGPSSATNUM; i++)
	{
		if (sat_rover[i].valid && sat_base[i].valid && i != refSatIndex_GPS)
		{
			//第一行和第二行
			for (int j = 0; j < doubleDiffNum_GPS; j++)
			{
				//一行伪距一行相位
				D(row, j) = Q_S_ref_P;
				D(row + doubleDiffNum_GPS, j + doubleDiffNum_GPS) = Q_S_ref_P;
				D(row + doubleDiffNum_GPS*2, j+ doubleDiffNum_GPS*2) = Q_S_ref_L;
				D(row + doubleDiffNum_GPS*3, j + doubleDiffNum_GPS*3) = Q_S_ref_L;
			}

			for (int j = 0; j < doubleDiffNum_GPS; j++)
			{
				if (row == j)
				{
					double Q_S_j_L = CalSigma2(obs_base.obsd[i].elev*DEG2RAD) +
						CalSigma2(obs_rover.obsd[i].elev*DEG2RAD);
					double Q_S_j_P = Q_S_j_L * multiplePL;
					D(row, j) = Q_S_ref_P + Q_S_j_P;
					D(row + doubleDiffNum_GPS, j + doubleDiffNum_GPS) = Q_S_ref_P + Q_S_j_P;
					D(row + doubleDiffNum_GPS * 2, j + doubleDiffNum_GPS * 2) = Q_S_ref_L + Q_S_j_L;
					D(row + doubleDiffNum_GPS * 3, j + doubleDiffNum_GPS * 3) = Q_S_ref_L + Q_S_j_L;
				}
			}		
			row++;
		}
	}


	//北斗方差阵构建
	int var_offset = 4 * obsNum_GPS - 4;	//BD方差阵的偏移量
	double Q_S_ref_L_BDS = CalSigma2(obs_base.obsd[refSatIndex_BDS].elev*DEG2RAD) +
		CalSigma2(obs_rover.obsd[refSatIndex_BDS].elev*DEG2RAD);
	double Q_S_ref_P_BDS = Q_S_ref_L_BDS * multiplePL;

	int doubleDiffNum_BDS = obsNum_BDS - 1;
	row = 0;
	for (int i = MAXGPSSATNUM; i < MAXSATNUM; i++)
	{
		if (sat_rover[i].valid && sat_base[i].valid && i != refSatIndex_BDS)
		{
			//第一行和第二行
			for (int j = 0; j < doubleDiffNum_BDS; j++)
			{
				//一行伪距一行相位
				D(var_offset+row, j + var_offset) = Q_S_ref_P_BDS;
				D(row + doubleDiffNum_BDS + var_offset, j + doubleDiffNum_BDS + var_offset) = Q_S_ref_P_BDS;
				D(row + doubleDiffNum_BDS * 2 + var_offset, j + doubleDiffNum_BDS * 2 + var_offset) = Q_S_ref_L_BDS;
				D(row + doubleDiffNum_BDS * 3 + var_offset, j + doubleDiffNum_BDS * 3 + var_offset) = Q_S_ref_L_BDS;
			}

			for (int j = 0; j < doubleDiffNum_BDS; j++)
			{
				if (row == j)
				{
					double Q_S_j_L_BDS = CalSigma2(obs_base.obsd[i].elev*DEG2RAD) +
						CalSigma2(obs_rover.obsd[i].elev*DEG2RAD);
					double Q_S_j_P_BDS = Q_S_j_L_BDS * multiplePL;
					D(row + var_offset, j + var_offset) = Q_S_ref_P_BDS + Q_S_j_P_BDS;
					D(row + doubleDiffNum_BDS + var_offset, j + doubleDiffNum_BDS + var_offset) = Q_S_ref_P_BDS + Q_S_j_P_BDS;
					D(row + doubleDiffNum_BDS * 2 + var_offset, j + doubleDiffNum_BDS * 2 + var_offset) = Q_S_ref_L_BDS + Q_S_j_L_BDS;
					D(row + doubleDiffNum_BDS * 3 + var_offset, j + doubleDiffNum_BDS * 3 + var_offset) = Q_S_ref_L_BDS + Q_S_j_L_BDS;
				}
			}
			row++;
		}
	}

	//B.Print();
	//D.Print();
	//l.Print();
	
	P = D.inverse();
	Matrix x(3 + 2 * obsNum - 4, 1);
	Matrix v(4 * obsNum - 8, 1);

	x = (B.transpose()*P*B).inverse()*B.transpose()*P*l;
	v = B*x - l;
	Matrix Q(3 + 2 * obsNum_GPS - 4, 3 + 2 * obsNum - 4);

	//x.Print();

	Q = (B.transpose()*P*B).inverse();	
	rtk_float.Q = Q;
	rtk_float.x = x;
	for (int i = 0; i < 3; i++)
	{
		rtk_float.pos[i] = sol_rover.pos[i] + x(i, 0);
		rtk_float._baseline_float[i] = rtk_float.pos[i] - kPreciseBaseCoor[i];
	}


	//将浮点解坐标改正数加到sol_
	for (int i = 0; i < 3; i++)
	{
		sol_rover.pos[i] += x(i, 0);
	}


	int sod = obs_base.ct.Hour * 3600 + obs_base.ct.Minute * 60 + static_cast<int>(obs_base.ct.Second);
	if (sod == 46910+1800) //输出一个小时 0406:30600-start 0514:46910
		throw runtime_error("bug");	
}

/***********************************************************
 * func    ambiguity fixing module
 * @paras  [in]  Q:浮点解最小二乘的方差阵
 *				 x:浮点解最小二乘的待估参数，包括xyz改正数和模糊度浮点解
 *				 rtk_float:存放rtk浮点解的结构体变量
 *         [out] p_file:文件指针，将结果输出至文件
 *				 rtk_fix:存放rtk固定解的结构体变量
 *				 
 * @return void
***********************************************************/
void FixAmbiguity(RtkFloat& rtk_float, RtkFix& rtk_fix)
{
	//为调用lambda函数生成必要的参数
	int num_ambiguity = rtk_float.x.GetRows() - 3;	//模糊度个数
	int m = 2;								//固定解的组数
	double* a = new double[num_ambiguity];	//模糊度浮点解
	double* Qa = new double[num_ambiguity*num_ambiguity];	//模糊度的方差
	double* F = new double[num_ambiguity*m];//模糊度固定解
	double* s = new double[m];				//sum of squared residuals of fixed solutions

	for (int i = 0; i < num_ambiguity; i++)
		a[i] = rtk_float.x(i + 3, 0);
	for (int i = 0; i < num_ambiguity; i++)
		for (int j = 0; j < num_ambiguity; j++)
			Qa[j*num_ambiguity + i] = rtk_float.Q(i + 3, j + 3);
	//调用lambda函数
	int res = lambda(num_ambiguity, 2, a, Qa, F, s);
	double ratio = s[1] / s[0];
	rtk_fix.ratio = ratio;
	//rtk_fix.Qaa = Matrix(num_ambiguity, num_ambiguity);
	//for (int i = 0; i < num_ambiguity; i++)
	//{
	//	for (int j = 0; j < num_ambiguity; j++)
	//	{
	//		rtk_fix.Qaa(i, j) = F[i + j*num_ambiguity];
	//	}
	//}

	//计算基线固定解
	Matrix Qba(3, num_ambiguity);
	Matrix Qab(num_ambiguity, 3);
	Matrix Qbb_float(3, 3);
	Matrix Qaa(num_ambiguity, num_ambiguity);
	Matrix delta_a(num_ambiguity, 1);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Qbb_float(i, j) = rtk_float.Q(i, j);
	for (int i = 0; i < num_ambiguity; i++)
		for (int j = 0; j < num_ambiguity; j++)
			Qaa(i, j) = rtk_float.Q(i + 3, j + 3);
	for (int i = 0; i < num_ambiguity; i++)
	{
		Qab(i, 0) = rtk_float.Q(i + 3, 0);
		Qab(i, 1) = rtk_float.Q(i + 3, 1);
		Qab(i, 2) = rtk_float.Q(i + 3, 2);
	}
	for (int i = 0; i < num_ambiguity; i++)
	{
		Qba(0, i) = rtk_float.Q(0, i + 3);
		Qba(1, i) = rtk_float.Q(1, i + 3);
		Qba(2, i) = rtk_float.Q(2, i + 3);
	}
	for (int i = 0; i < num_ambiguity; i++)
		delta_a(i, 0) = rtk_float.x(i + 3, 0) - F[i];
	Matrix b_fix(3, 1);
	Matrix b_float(3, 1);
	for (int i = 0; i < 3; i++)
	{
		b_float(i, 0) = rtk_float._baseline_float[i];
	}
	//得到固定解 及其 方差阵
	b_fix = b_float - Qba*Qaa.inverse()*delta_a;
	Matrix Qbb_fix(3, 3);
	Qbb_fix = Qbb_float - Qba*Qaa.inverse()*Qab;
	for (int i = 0; i < 3; i++)
	{
		rtk_fix._baseline_fix[i] = b_fix(i, 0);
	}

	delete[] a;
	delete[] Qa;
	delete[] F;
	delete[] s;

}


//高度角法算sigma^2, E是高度角，单位rad
double CalSigma2(double E)
{
	double a = 0.004;
	double b = 0.003;
	return SQR(a) + SQR(b) / SQR(sin(E));
}