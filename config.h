#pragma once
#include <string>

extern const std::string kEphFilepath;	//星历的文件路径
extern const std::string kBaseFilepath;	//基准站数据的文件路径
extern const std::string kRoverFilepath;//流动站数据的文件路径
extern const char kFloatResOutFilepath[60];//浮点解输出的文件路径
extern const char kFixResOutFilepath[60];	//固定解输出的文件路径
extern const double kPreciseRoverCoor[3];	//基准站的精确坐标
extern const double kPreciseBaseCoor[3];	//流动站的精确坐标

extern const double kElevationMaskAngle;	//截止高度角，单位:deg
