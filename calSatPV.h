#pragma once
#include "global.h"
#include "rinex.h"



//����������tʱ�̵�PV
void calSingleSatPV(MJD imputT, const eph_t& eph, sat_t& sat);
void CalBatchSatPV(MJD inputT, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM]);
void calSatTransPV(obs_t& obst, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM]);
void calSatTransPVBDS(obs_t& obst, eph_t eph[MAXSATNUM][EPOCHNUM], sat_t sat[MAXSATNUM]);