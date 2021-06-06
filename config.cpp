#include "config.h"

const std::string kEphFilepath = "NovAtel短基线\\0514\\brdm1340.21p";	//星历的文件路径
const std::string kBaseFilepath = "NovAtel短基线\\0514\\SGG_20210514.21O";	//基准站数据的文件路径
const std::string kRoverFilepath = "NovAtel短基线\\0514\\CENT_20210514.21O";//流动站数据的文件路径
const char kFloatResOutFilepath[60] = "Resutls\\NovATel短基线0514\\resfloat.txt";//浮点解输出的文件路径
const char kFixResOutFilepath[60] = "Resutls\\NovATel短基线0514\\resfix.txt";	//固定解输出的文件路径

//0514 short baseline
const double kPreciseRoverCoor[3] = { -2267810.910,5009435.734,3220956.013 };	//基准站的精确坐标
const double kPreciseBaseCoor[3] = { -2267810.196,5009356.572,3221000.818 };	//流动站的精确坐标

//0406 zero baseline
//const double kPreciseBaseCoor[3] = {-2267810.196,5009356.572,3221000.818 };
//const double kPreciseRoverCoor[3] = {-2267810.196,5009356.572,3221000.818 };

const double kElevationMaskAngle = 15.0;	//截止高度角，单位：deg