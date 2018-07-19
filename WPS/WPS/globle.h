
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include <math.h>
#include "..\..\BLPK\blp.h"

#define STR_LEN 32
#define POINTS_PER_DAY 24
#define COLS (POINTS_PER_DAY+5)
typedef char str[STR_LEN+1];
#define STR_LEN_LONG 1024
typedef char str_long[STR_LEN_LONG+1];
#define DAYS_PER_YEAR 365
#define WPOS_MAX_PROJECT 100 /*�������󷽰���*/

/*pa��һ��ͨ��Ŀ¼�Ľṹ��n��ʾĿ¼�к��е���Ŀ����cap��ʾĿ¼�ܹ����ɵ���Ŀ����
** buf��ָ��ָ�������ָ�룬buf[i]��ʾĿ¼�е�i����Ŀ�ĵ�ַ����һ�㶼ָ������ϵͳԪ�������ԵĽṹ��*/
typedef struct pa
{
	int n;
	int cap;
	void** buf;
}pa;

/*�������*/
typedef struct unit_para
{
	int nDate,nDate2;/*��ʼ����*/
	int nNodeId;/*�����ڵ�ID*/
	int nNodeIdx;/*�����ڵ㣭�����ڲ����*/
	int nPlantId;/*�����糧ID*/
	int nCompanyID;/*�������繫˾ID*/
	int nUnitTypeIn;/*����Ļ�������,1:���;2:��ͣ;3:E��ȼ��;4:�˵�;5:ˮ��;6:����;7:����; 8�ȵ磻9:F��ȼ����10:110kVͳ����11:��磻20:��ͳ��zhangning20121207*/
	double dPowerMax;/*������*/
	double dPowerMin;/*��С����*/

	double dUpDownCost;/*�������ͣ����*/
	double dEnergyCapa;/*ˮ����������������ˮ��*/
	double dConvertEfft;/*�����ת��Ч��*/
	double dFOR;/*�����ǿ��ͣ����*/
	double dSelfUse;/*���õ���*/
	double dPowerEffect;/*���Ϊ����ฺ�ɣ���Ϊ1�����Ϊ���ɲฺ�ɣ���Ϊ=1-���õ���-����*/
	double dVariableCost;/*��߳�����Ӧ�Ŀɱ����гɱ�*/
	double dVariableCost2;/*��ͳ�����Ӧ�Ŀɱ����гɱ�*/
	double dCoal;/*��߳�����Ӧ��ú��*/
	double dCoal2;/*��߳�����Ӧ��ú��*/
	double dRampRate;/*�����ٶȣ�MW/min��*/
}unit_para;
/*��糡����*/
typedef struct wind_para
{
	int nDate,nDate2;/*��ʼ����*/
	int nNumofTurbine;//�������
	double dCutinSP;//�������
	double dRatedSP;//�����
	double dCutoutSP;//�г�����
	double dWakeeffect;//β��ϵ��
	double dAvailableRate;//������
	double dForecastError;//Ԥ��������ռװ�������İٷֱ�
	int nWindZone;//���ڷ���Id
}wind_para;


/*����ָ������*/
typedef struct unit_power
{
	int nDate,nDate2;
	int nFixedPower;/*�������ճ���Ϊ�̶�����*/
	double dPower[POINTS_PER_DAY];
}unit_power;

/*����*/
typedef struct unit
{
	int nId;/*����ID*/
	int nIdx;/*�����������*/
	pa para;/*unit_para*/
	pa windpara;/*wind_para*/
	pa power;/*����ָ������*/
	pa fpower;/*����Ԥ�����*/
	pa appStateDate;/*ָ��״̬����*/
	pa appStateDate2;/*ָ��״̬��ֹ����*/
	pa appState;/*ָ��״̬*/
	pa mt1;/*������ʼ����*/
	pa mt2;/*������ֹ���ڣ����������գ�*/

#define UT_WIND 10/*������*/

	int nUnitType;
	int nStateOut;/*�������״̬*/
	int nStateIn;/*��������״̬*/
}unit;

typedef struct wind_zone
{
	int nId;/*����ID*/
	pa para;
}wind_zone;

typedef struct wind_zone_para
{
	int nDate,nDate2;
	double dWeibull_c;//����Weibull�ֲ�����c
	double dWeibull_k;//����Weibull�ֲ�����k
	double dtheta;//��������غ���˥��ϵ��
	double nDayAve[POINTS_PER_DAY];//��ƽ������
	double nMonthAve[12];//��ƽ������
}wind_zone_para;

typedef struct wind_powerPDF
{
	int nDate,nDate2;
	mtx *PowerPDF;
}wind_powerPDF;

typedef struct load
{
	int nDate,nDate2;
	double dLoad[POINTS_PER_DAY];
}load;

typedef struct regulation
{
	int nDate,nDate2;
	double dNegRegulationPro[12];
	int nNegRegulationdays[12];
	double dFlatRegulationPro[12];
	int nFlatRegulationdays[12];
	double dPosRegulationPro[12];
	int nPosRegulationdays[12];
}regulation;

typedef struct app
{
	int nPJID;
	int nStartDate,nEndDate;
	int nTimes;/*Ŀǰģ��ĳ�����*/
	double dWindCCFor;
	pa aUnits;/*unit*/
	//pa aNodes;/*node*/
	clock_t t_start,t_finish;

	pa aSys_load;  /*load*/ //�洢ϵͳ��������
	pa aRegulationPro;   /*regulation*/ //�洢�����������巴�������

	//ϵͳ����
	int bWindFcst; /*�Ƿ�ģ����Ԥ���������ģ������Ļ�����ģ��Ԥ�����*/
	int bCorrlation;/*�Ƿ��Ƿ�糡�ռ������*/
	int bReliability;/*�Ƿ��Ƿ�����Ŀɿ���*/
	int bWakeEffect;/*�Ƿ��Ƿ�糡β��ЧӦ*/
	int bDayRhythm;/*�Ƿ��Ƿ���������*/
	int bSeasonRhythm;/*�Ƿ��Ƿ��ټ�������*/
	int nWindPowerCurve;/*���������������ģʽ��0��ֱ���ͣ�1������������*/
	double dSpeedBias;/*ƽ�����ٵ�����*/
	int nErrorType; /*Ԥ���������ͣ�0������Ԥ����������̬�ֲ���1������Ԥ����������̬�ֲ�*/
	int nErrorCoherence; /*Ԥ�����һ���ԣ�0 ��ʱ��Ԥ������׼����ͬ��1����ʱ��Ԥ����������ӣ�������Ԥ�����ӵ�һ��12��Ԥ�⵽�ڶ�������24�㣬��ֱ��ģ�Ϳ���*/
	int nRandomType; /*������������ɷ�ʽ��0, ������ӣ�1���̶�����*/
	int nRandomSeed; /*���������*/
	int nRandomSeednow; /*��������Ӽ�¼�����ڶ�β���*/
	int nSimulateTimes; /*ģ�ⳡ����*/
	int bAppPower;//�Ƿ����ⲿ����ķ��ָ�����������ѡ0���򲻿���ָ���������г�����ģ�⣬���ѡ1������ָ��������ָ�������������ָ��������Ӧ��ʱ�ν���ģ��zhangning20130929
	double dPeakValleyRegulationLimit; //�����������巴�������ʱ�����ı�ϵͳ��Ȳ�İٷֱ����ơ������ðٷֱȾ���Ϊ��������򷴵��塣������Ϊ�Է�Ȳ�ûӰ�졣

	/*ÿ��ļ����־����ȫ�ֱ�־��ÿ����һ�춼���ܸı�*/
	//int nCUnits,nFUnits,nWUnits,nHUnits,nLines,nSections,nNodes; 
	int nWUnits; 
	int nCalcFlag;
	pa aWindZone;/*wind_zone*/
	mtx *pmWindCorr;/*�������ϵ������*/
	pa aWindPowerPDF;/*�������ֲ�*/

	FILE* fRet;

	char caExePathName[255];			/*RECC.exe�ļ����ļ�·��*/
	char caInPathName[255];             /*���������ļ����ļ�·��*/
	char caOutPathName[255];            /*��������ļ����ļ�·��*/ 

	char sDir[WPOS_MAX_PROJECT][255]; /*��Ŷ���Ŀ¼*/
	char sCursDir[255]; /*��ǰ����Ŀ¼*/
	int	 nDir; /*����Ŀ¼����*/

}app;

