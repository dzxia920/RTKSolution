#pragma once

#include <string>
#include <fstream>
#include <cstdlib>

#include "global.h"
#include "rinex.h"
#include "TimeCoor.h"

extern int OBSTYPESEQUENCE[MAXOBSSYSTEM][8];

// ��ȡ�۲�ֵ
void readObs(std::fstream& fs, obs_t &obs);
void readObsHeader(std::fstream& fs, obs_header &obsHeader);

// ��ȡ��������
void readNav(std::fstream& fin,eph_t eph[MAXSATNUM][EPOCHNUM]);
void readNavHeader(std::fstream& fs, nav_header &navHeader);