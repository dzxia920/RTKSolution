#include "config.h"

const std::string kEphFilepath = "NovAtel�̻���\\0514\\brdm1340.21p";	//�������ļ�·��
const std::string kBaseFilepath = "NovAtel�̻���\\0514\\SGG_20210514.21O";	//��׼վ���ݵ��ļ�·��
const std::string kRoverFilepath = "NovAtel�̻���\\0514\\CENT_20210514.21O";//����վ���ݵ��ļ�·��
const char kFloatResOutFilepath[60] = "Resutls\\NovATel�̻���0514\\resfloat.txt";//�����������ļ�·��
const char kFixResOutFilepath[60] = "Resutls\\NovATel�̻���0514\\resfix.txt";	//�̶���������ļ�·��

//0514 short baseline
const double kPreciseRoverCoor[3] = { -2267810.910,5009435.734,3220956.013 };	//��׼վ�ľ�ȷ����
const double kPreciseBaseCoor[3] = { -2267810.196,5009356.572,3221000.818 };	//����վ�ľ�ȷ����

//0406 zero baseline
//const double kPreciseBaseCoor[3] = {-2267810.196,5009356.572,3221000.818 };
//const double kPreciseRoverCoor[3] = {-2267810.196,5009356.572,3221000.818 };

const double kElevationMaskAngle = 15.0;	//��ֹ�߶Ƚǣ���λ��deg