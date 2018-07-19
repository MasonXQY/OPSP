
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
#define WPOS_MAX_PROJECT 100 /*计算的最大方案数*/

/*pa是一个通用目录的结构，n表示目录中含有的条目数，cap表示目录能够容纳的条目数，
** buf是指向指针数组的指针，buf[i]表示目录中第i个条目的地址，他一般都指向描述系统元件的属性的结构体*/
typedef struct pa
{
	int n;
	int cap;
	void** buf;
}pa;

/*机组参数*/
typedef struct unit_para
{
	int nDate,nDate2;/*起始日期*/
	int nNodeId;/*所处节点ID*/
	int nNodeIdx;/*所处节点－程序内部编号*/
	int nPlantId;/*所属电厂ID*/
	int nCompanyID;/*所属发电公司ID*/
	int nUnitTypeIn;/*输入的机组类型,1:火电;2:启停;3:E级燃机;4:核电;5:水电;6:抽蓄;7:区外; 8热电；9:F级燃机；10:110kV统调；11:风电；20:非统调zhangning20121207*/
	double dPowerMax;/*最大出力*/
	double dPowerMin;/*最小出力*/

	double dUpDownCost;/*机组的启停费用*/
	double dEnergyCapa;/*水电容量或抽蓄的最大抽水量*/
	double dConvertEfft;/*抽蓄的转换效率*/
	double dFOR;/*机组的强迫停运率*/
	double dSelfUse;/*厂用电率*/
	double dPowerEffect;/*如果为发电侧负荷，则为1，如果为负荷侧负荷，则为=1-厂用电率-线损*/
	double dVariableCost;/*最高出力对应的可变运行成本*/
	double dVariableCost2;/*最低出力对应的可变运行成本*/
	double dCoal;/*最高出力对应的煤耗*/
	double dCoal2;/*最高出力对应的煤耗*/
	double dRampRate;/*爬坡速度（MW/min）*/
}unit_para;
/*风电场参数*/
typedef struct wind_para
{
	int nDate,nDate2;/*起始日期*/
	int nNumofTurbine;//风机数量
	double dCutinSP;//切入风速
	double dRatedSP;//额定风速
	double dCutoutSP;//切出风速
	double dWakeeffect;//尾流系数
	double dAvailableRate;//可用率
	double dForecastError;//预测绝对误差占装机容量的百分比
	int nWindZone;//所在风区Id
}wind_para;


/*机组指定出力*/
typedef struct unit_power
{
	int nDate,nDate2;
	int nFixedPower;/*标明该日出力为固定出力*/
	double dPower[POINTS_PER_DAY];
}unit_power;

/*机组*/
typedef struct unit
{
	int nId;/*机组ID*/
	int nIdx;/*机组索引编号*/
	pa para;/*unit_para*/
	pa windpara;/*wind_para*/
	pa power;/*机组指定出力*/
	pa fpower;/*机组预测出力*/
	pa appStateDate;/*指定状态日期*/
	pa appStateDate2;/*指定状态终止日期*/
	pa appState;/*指定状态*/
	pa mt1;/*检修起始日期*/
	pa mt2;/*检修终止日期（不包含本日）*/

#define UT_WIND 10/*风电机组*/

	int nUnitType;
	int nStateOut;/*机组输出状态*/
	int nStateIn;/*机组输入状态*/
}unit;

typedef struct wind_zone
{
	int nId;/*风区ID*/
	pa para;
}wind_zone;

typedef struct wind_zone_para
{
	int nDate,nDate2;
	double dWeibull_c;//风速Weibull分布参数c
	double dWeibull_k;//风速Weibull分布参数k
	double dtheta;//风速自相关函数衰减系数
	double nDayAve[POINTS_PER_DAY];//日平均风速
	double nMonthAve[12];//月平均风速
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
	int nTimes;/*目前模拟的场景数*/
	double dWindCCFor;
	pa aUnits;/*unit*/
	//pa aNodes;/*node*/
	clock_t t_start,t_finish;

	pa aSys_load;  /*load*/ //存储系统负荷曲线
	pa aRegulationPro;   /*regulation*/ //存储风电各月正调峰反调峰概率

	//系统参数
	int bWindFcst; /*是否模拟风电预测出力，在模拟出力的基础上模拟预测出力*/
	int bCorrlation;/*是否考虑风电场空间相关性*/
	int bReliability;/*是否考虑风电机组的可靠性*/
	int bWakeEffect;/*是否考虑风电场尾流效应*/
	int bDayRhythm;/*是否考虑风速日特性*/
	int bSeasonRhythm;/*是否考虑风速季节特性*/
	int nWindPowerCurve;/*风机功率特性曲线模式，0，直线型，1，三次曲线型*/
	double dSpeedBias;/*平均风速调整量*/
	int nErrorType; /*预测误差的类型，0，出力预测误差服从正态分布，1，风速预测误差服从正态分布*/
	int nErrorCoherence; /*预测误差一致性：0 各时段预测误差标准差相同，1，各时段预测误差逐渐增加，这里面预测误差按从第一天12点预测到第二天晚上24点，按直线模型考虑*/
	int nRandomType; /*随机数种子生成方式：0, 随机种子，1，固定种子*/
	int nRandomSeed; /*随机数种子*/
	int nRandomSeednow; /*随机数种子记录，用于多次采样*/
	int nSimulateTimes; /*模拟场景数*/
	int bAppPower;//是否考虑外部输入的风电指定出力，如果选0，则不考虑指定出力所有出力均模拟，如果选1，则有指定出力按指定出力输出，无指定出力对应的时段进行模拟zhangning20130929
	double dPeakValleyRegulationLimit; //计算风电正调峰反调峰概率时，风电改变系统峰谷差的百分比限制。超出该百分比就认为是正调峰或反调峰。其余认为对峰谷差没影响。

	/*每天的计算标志，非全局标志，每计算一天都可能改变*/
	//int nCUnits,nFUnits,nWUnits,nHUnits,nLines,nSections,nNodes; 
	int nWUnits; 
	int nCalcFlag;
	pa aWindZone;/*wind_zone*/
	mtx *pmWindCorr;/*风速相关系数矩阵*/
	pa aWindPowerPDF;/*风电出力分布*/

	FILE* fRet;

	char caExePathName[255];			/*RECC.exe文件的文件路径*/
	char caInPathName[255];             /*输入数据文件的文件路径*/
	char caOutPathName[255];            /*输出数据文件的文件路径*/ 

	char sDir[WPOS_MAX_PROJECT][255]; /*存放读入目录*/
	char sCursDir[255]; /*当前运行目录*/
	int	 nDir; /*读入目录个数*/

}app;

