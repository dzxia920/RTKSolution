#pragma once
#include "global.h"
#include "rinex.h"

/*
*	�������˵����
*	���ڼ����źŷ���ʱ����λ���Լ������������ЧӦ���ʴ˴����ٸ�
*	�������ӳ٣�Hopfieldģ��
*	������ӳ٣�������������˫Ƶα�����
*	������ת����������P129 �൱������������ת��һ���Ƕ�alpha��
*	TGD����û��
*/

void correctEarthRot(sat_t& sat);


//�����߶Ƚǡ㣬��վ�̣߳�����Hopfieldģ�ͼ���������ӳٸ���
double Hopfield(double E, double hs);
