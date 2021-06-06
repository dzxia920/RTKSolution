#include "readRinex.h"

using namespace std;

//*************函数内部使用的常/变量****************
const string EMPTYP(14, ' ');
int OBSTYPESEQUENCE[MAXOBSSYSTEM][8] = { 0 };
static int CurrentEPOCHNUM[MAXSATNUM] = { 0 };//指示EPOCHNUM现在有多少历元的数组

//*************函数内部使用的小函数*****************
//可以考虑用scanf实现
MJD timeStr2MJD(const string& timeStamp);
CommonTime timeStr2CommonTime(const string& timeStamp);
MJD dateStr2MJD(const std::string& date);


void readObsHeader(std::fstream& fs, obs_header &obsHeader)
{
	string headLine;		//存放卫星数据的行
	
	//读取header
	do {
		getline(fs, headLine);
		//如果有接收机的概略位置，则读取
		if (headLine.substr(60, 19) == "APPROX POSITION XYZ")
		{
			obsHeader.approxpos_rcv_xyz[0] = stod(headLine.substr(0, 14));
			obsHeader.approxpos_rcv_xyz[1] = stod(headLine.substr(14, 14));
			obsHeader.approxpos_rcv_xyz[2] = stod(headLine.substr(28, 14));
		}

		//读取采样间隔
		if (headLine.substr(60, 8) == "INTERVAL")
			obsHeader.interval = stod(headLine.substr(5, 5));
		

		//读取SYS/#/OBS TYPES，将其转成顺序
		int obsTypeNum = 0;	//每个系统观测值类型的数目 e.g.C1C
		while (headLine.substr(60, 19) == "SYS / # / OBS TYPES")
		{
			//三系统以上时可以更完善
			switch (headLine.at(0))
			{
			case 'G':
				obsTypeNum = stoi(headLine.substr(3, 3));
				obsHeader.obsType[0] = headLine.substr(0, 60);
				if (obsTypeNum > 13)
				{
					getline(fs, headLine);
					obsHeader.obsType[0] += headLine.substr(7, 53);
				}
				break;
			case 'C':
				obsTypeNum = stoi(headLine.substr(3, 3));
				obsHeader.obsType[1] = headLine.substr(0, 60);
				if (obsTypeNum > 13)
				{
					getline(fs, headLine);
					obsHeader.obsType[1] += headLine.substr(7, 53);
				}
				break;
			default:
				break;
			}
			getline(fs, headLine);
		}
		obsHeader.calobsTypeSequence(OBSTYPESEQUENCE);
	} while (headLine.substr(60, 13) != "END OF HEADER");
}

void readObs(std::fstream& fs, obs_t &obst)
{
	string satline;		//存放卫星数据的行
	string timeStamp;	//每个历元的时间戳行

	int SatNumofCurrentEpoch = 0;		//当前历元可见卫星数

	int SatNumofGC = 0;					//GPS BDS Sat Num

	//将所有卫星的valid标志置成false
	for (int i = 0; i < MAXSATNUM; i++)
	{
		obst.obsd[i].valid = false;
	}

	if (fs.peek() == '>')
	{
		getline(fs, timeStamp);
		//tline 转成时间 MJD？ 读取卫星行数
		SatNumofCurrentEpoch = stoi(timeStamp.substr(32, 3));
		obst.mjd = timeStr2MJD(timeStamp);
		obst.ct = timeStr2CommonTime(timeStamp);

		//for loop 读取这么多行的卫星数据 
		// 根据RPN号存储，0-31GPS  32-93BDS
		int prnIndex_GPS = 0;
		int prnIndex_BDS = 0;
		for (int i = 0; i < SatNumofCurrentEpoch; i++)
		{
			getline(fs, satline);
			char sys = satline[0];
			switch (sys)
			{
			case 'G': 
				if (satline.substr(3 + OBSTYPESEQUENCE[0][0] * 16, 14) == EMPTYP
					|| satline.substr(3 + OBSTYPESEQUENCE[0][1] * 16, 14) == EMPTYP
					|| satline.substr(3 + OBSTYPESEQUENCE[0][4] * 16, 14) == EMPTYP
					|| satline.substr(3 + OBSTYPESEQUENCE[0][5] * 16, 14) == EMPTYP)
					break;
				prnIndex_GPS = stoi(satline.substr(1, 2)) - 1;	//PRN -1
				obst.obsd[prnIndex_GPS].valid = true;
				obst.obsd[prnIndex_GPS].sys = sys;
				obst.obsd[prnIndex_GPS].prn = stoi(satline.substr(1, 2));
				obst.obsd[prnIndex_GPS].P[0] = stod(satline.substr(3 + OBSTYPESEQUENCE[0][0] * 16, 14));
				obst.obsd[prnIndex_GPS].L[0] = stod(satline.substr(3 + OBSTYPESEQUENCE[0][1] * 16, 14));
				obst.obsd[prnIndex_GPS].D[0] = stod(satline.substr(3 + OBSTYPESEQUENCE[0][2] * 16, 14));
				obst.obsd[prnIndex_GPS].SNR[0] = stod(satline.substr(3 + OBSTYPESEQUENCE[0][3] * 16, 14));
				obst.obsd[prnIndex_GPS].P[1] = stod(satline.substr(3 + OBSTYPESEQUENCE[0][4] * 16, 14));
				obst.obsd[prnIndex_GPS].L[1] = stod(satline.substr(3 + OBSTYPESEQUENCE[0][5] * 16, 14));
				obst.obsd[prnIndex_GPS].D[1] = stod(satline.substr(3 + OBSTYPESEQUENCE[0][6] * 16, 14));
				obst.obsd[prnIndex_GPS].SNR[1] = stod(satline.substr(3 + OBSTYPESEQUENCE[0][7] * 16, 14));
				SatNumofGC++;
				break;
			case 'C':
				if (satline.substr(3 + OBSTYPESEQUENCE[1][0] * 16, 14) == EMPTYP
					|| satline.substr(3 + OBSTYPESEQUENCE[1][1] * 16, 14) == EMPTYP
					|| satline.substr(3 + OBSTYPESEQUENCE[1][4] * 16, 14) == EMPTYP
					|| satline.substr(3 + OBSTYPESEQUENCE[1][5] * 16, 14) == EMPTYP)
					break;
				prnIndex_BDS = stoi(satline.substr(1, 2)) - 1 + 32;
				////不要GEO卫星
				//if (prnIndex_BDS + 1 - 32 < 6 || prnIndex_BDS + 1 - 32 > 58)
				//	break;
				obst.obsd[prnIndex_BDS].valid = true;
				obst.obsd[prnIndex_BDS].sys = sys;
				obst.obsd[prnIndex_BDS].prn = stoi(satline.substr(1, 2));
				obst.obsd[prnIndex_BDS].P[0] = stod(satline.substr(3 + OBSTYPESEQUENCE[1][0] * 16, 14));
				obst.obsd[prnIndex_BDS].L[0] = stod(satline.substr(3 + OBSTYPESEQUENCE[1][1] * 16, 14));
				obst.obsd[prnIndex_BDS].D[0] = stod(satline.substr(3 + OBSTYPESEQUENCE[1][2] * 16, 14));
				obst.obsd[prnIndex_BDS].SNR[0] = stod(satline.substr(3 + OBSTYPESEQUENCE[1][3] * 16, 14));
				obst.obsd[prnIndex_BDS].P[1] = stod(satline.substr(3 + OBSTYPESEQUENCE[1][4] * 16, 14));
				obst.obsd[prnIndex_BDS].L[1] = stod(satline.substr(3 + OBSTYPESEQUENCE[1][5] * 16, 14));
				obst.obsd[prnIndex_BDS].D[1] = stod(satline.substr(3 + OBSTYPESEQUENCE[1][6] * 16, 14));
				obst.obsd[prnIndex_BDS].SNR[1] = stod(satline.substr(3 + OBSTYPESEQUENCE[1][7] * 16, 14));
				SatNumofGC++;
				break;
			default:
				break;
			}
		}

	}
	else getline(fs, satline);
	obst.n = SatNumofGC;
}

void readNavHeader(std::fstream& fin, nav_header &navHeader)
{
	string headerLine;
	//文件流的状态，用以指示是否读到文件头末尾
	bool fileState = true;

	for (size_t i = 0; getline(fin, headerLine) && fileState; i++)
	{
		if (headerLine.substr(60, 13) == "END OF HEADER")
		{
			fileState = false;
			break;
		}
	}
}

//NEED to IMPROVE
void readNav(std::fstream& fin,eph_t eph[MAXSATNUM][EPOCHNUM])
{
	string oribit[8];
	int i = 0;
	do {		
		i = 0;		
		//根据星历文件格式读取参数，看系统是否匹配
		int PRNIndex_GPS = 0, PRNIndex_BDS = 0;
		
		switch (fin.peek())
		{
		case 'C':
			while (i != 8 && getline(fin, oribit[i]))
				i++;
			PRNIndex_BDS = stoi(oribit[0].substr(1, 2)) + MAXGPSSATNUM - 1;
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].valid = true;
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].sys = oribit[0].at(0);
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].prn = stoi(oribit[0].substr(1, 2));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].toc = dateStr2MJD(oribit[0].substr(4, 19));	
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].a0 = stod(oribit[0].substr(23, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].a1 = stod(oribit[0].substr(42, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].a2 = stod(oribit[0].substr(61, 19));
											  
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].crs = stod(oribit[1].substr(23, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].delta_n = stod(oribit[1].substr(42, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].M0 = stod(oribit[1].substr(61, 19));
											  
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].cuc = stod(oribit[2].substr(4, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].e = stod(oribit[2].substr(23, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].cus = stod(oribit[2].substr(42, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].sqrtA = stod(oribit[2].substr(61, 19));

			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].sow = stod(oribit[3].substr(4, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].cic = stod(oribit[3].substr(23, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].omega0 = stod(oribit[3].substr(42, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].cis = stod(oribit[3].substr(61, 19));
				
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].i0 = stod(oribit[4].substr(4, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].crc = stod(oribit[4].substr(23, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].omega = stod(oribit[4].substr(42, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].omegaDot = stod(oribit[4].substr(61, 19));
				
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].idot = stod(oribit[5].substr(4, 19));
			eph[PRNIndex_BDS][CurrentEPOCHNUM[PRNIndex_BDS]].week = stoi(oribit[5].substr(42, 19));
			CurrentEPOCHNUM[PRNIndex_BDS]++;
			break;
		case 'G':
			while (i != 8 && getline(fin, oribit[i]))
				i++;
			PRNIndex_GPS = stoi(oribit[0].substr(1, 2)) - 1;
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].valid = true;
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].sys = oribit[0].at(0);
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].prn = stoi(oribit[0].substr(1, 2));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].toc = dateStr2MJD(oribit[0].substr(4, 19));
				
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].a0 = stod(oribit[0].substr(23, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].a1 = stod(oribit[0].substr(42, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].a2 = stod(oribit[0].substr(61, 19));
				
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].crs = stod(oribit[1].substr(23, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].delta_n = stod(oribit[1].substr(42, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].M0 = stod(oribit[1].substr(61, 19));
				
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].cuc = stod(oribit[2].substr(4, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].e = stod(oribit[2].substr(23, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].cus = stod(oribit[2].substr(42, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].sqrtA = stod(oribit[2].substr(61, 19));
				
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].sow = stod(oribit[3].substr(4, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].cic = stod(oribit[3].substr(23, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].omega0 = stod(oribit[3].substr(42, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].cis = stod(oribit[3].substr(61, 19));
				
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].i0 = stod(oribit[4].substr(4, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].crc = stod(oribit[4].substr(23, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].omega = stod(oribit[4].substr(42, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].omegaDot = stod(oribit[4].substr(61, 19));
				
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].idot = stod(oribit[5].substr(4, 19));
			eph[PRNIndex_GPS][CurrentEPOCHNUM[PRNIndex_GPS]].week = stoi(oribit[5].substr(42, 19));
			CurrentEPOCHNUM[PRNIndex_GPS]++;
			break;
		case 'E':
			while (i != 8 && getline(fin, oribit[i]))
				i++;
			break;
		case 'R':
			while (i != 4 && getline(fin, oribit[i]))
				i++;
			break;
		case 'S':
			while (i != 4 && getline(fin, oribit[i]))
				i++;
			break;
		default:
			getline(fin, oribit[0]);
			break;
		}		
	} while (!fin.eof());
}


MJD timeStr2MJD(const string& timeStamp)
{
	CommonTime ct;
	ct.Year = stoi(timeStamp.substr(2, 4));
	ct.Month = stoi(timeStamp.substr(7, 2));
	ct.Day = stoi(timeStamp.substr(10, 2));
	ct.Hour = stoi(timeStamp.substr(13, 2));
	ct.Minute = stoi(timeStamp.substr(16, 2));
	ct.Second = stod(timeStamp.substr(19, 11));

	MJD mjd = ct.ToMJD();
	return mjd;
}

CommonTime timeStr2CommonTime(const string& timeStamp)
{
	CommonTime ct;
	ct.Year = stoi(timeStamp.substr(2, 4));
	ct.Month = stoi(timeStamp.substr(7, 2));
	ct.Day = stoi(timeStamp.substr(10, 2));
	ct.Hour = stoi(timeStamp.substr(13, 2));
	ct.Minute = stoi(timeStamp.substr(16, 2));
	ct.Second = stod(timeStamp.substr(19, 11));
	return ct;
}

MJD dateStr2MJD(const std::string& date)
{
	CommonTime ct;
	MJD mjd;
	ct.Year = stoi(date.substr(0, 4));
	ct.Month = stoi(date.substr(5, 2));
	ct.Day = stoi(date.substr(8, 2));
	ct.Hour = stoi(date.substr(11, 2));
	ct.Minute = stoi(date.substr(14, 2));
	ct.Second = stod(date.substr(17, 2));
	mjd = ct.ToMJD();
	return mjd;
}

