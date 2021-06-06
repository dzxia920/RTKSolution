#pragma once
#include <string>
//全局变量头文件

#define MAXSATNUM 93
#define EPOCHNUM 30		//导航电文中每颗卫星可存储的最大星历数
#define MAXOBSSYSTEM 2	//最多的系统数，用于obs Header
#define MAXGPSSATNUM 32

#define PI 3.1415926535897932
#define CLIGHT 299792458.0
#define OMGE 7.2921151467E-5   /* earth angular velocity (IS-GPS) (rad/s) */
#define mu 3.986005E+14			//mu=GM,地球引力常量

#define MU_BDS 3.986004418E+14	//北斗参考椭球的常量
#define OMGE_BDS 7.2921150E-5
#define PI_BDS 3.1415926535898

#define GPSL1FREQUENCE 1575.42E6
#define GPSL2FREQUENCE 1227.60E6
#define BDSB1FREQUENCE 1561.098E6
#define BDSB3FREQUENCE 1268.052E6

#define GPSL1LAMBDA (CLIGHT/GPSL1FREQUENCE)
#define GPSL2LAMBDA (CLIGHT/GPSL2FREQUENCE)
#define BDSB1LAMBDA (CLIGHT/BDSB1FREQUENCE)
#define BDSB3LAMBDA (CLIGHT/BDSB3FREQUENCE)

#define DEG2RAD (PI/180.0)
#define RAD2DEG (180.0/PI)

#define SQR(x) ((x)*(x))  //平方

extern const std::string GPSOBSTYPE[8];
extern const std::string BDSOBSTYPE[8];

