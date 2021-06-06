#pragma once
#include "config.h"
#include "MatrixD.h"
#include "rinex.h"
#include "SPPSolution.h"
#include "lambda.h"
//建立函数模型
void CreateRTKModel(obs_t& obs_base, obs_t& obs_rover, sat_t sat_base[MAXSATNUM], sat_t sat_rover[MAXSATNUM], sol_t& sol_base, sol_t& sol_rover, RtkFloat& rtk_float);

//fix ambiguity
void FixAmbiguity(RtkFloat& rtk_float, RtkFix& rtk_fix);