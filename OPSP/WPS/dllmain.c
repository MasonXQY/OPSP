// dllmain.cpp : 定义 DLL 应用程序的入口点。
#include "stdafx.h"
#include "wps.h"
#include <stdio.h>
#include <stdlib.h>
#include <direct.h>                                                                                                        
#include <time.h>
#include <string.h>
#include <stdarg.h>
#include "..\BLPK\blp.h"
#include "globle.h"
#include "func.h"
#include "gConstant.h"

#ifdef WIN32
#include <Windows.h>
#endif
#pragma warning(disable: 4996)

#if _MSC_VER>=1900
#include "stdio.h" 
_ACRTIMP_ALT FILE* __cdecl __acrt_iob_func(unsigned);
#ifdef __cplusplus 
extern "C"
#endif 
FILE* __cdecl __iob_func(unsigned i) {
	return __acrt_iob_func(i);
}
#endif /* _MSC_VER>=1900 */


pa* pa_create()
{
	pa *p = malloc(sizeof(pa));
	memset(p, 0, sizeof(pa));
	return p;
}
void pa_clear(pa *p, int bFreeItem)
{
	int i;
	if (p == NULL)return;
	if (bFreeItem != 0)
	{
		for (i = 0; i<p->n; i++)free(p->buf[i]);
	}
	if (p->buf != NULL)
	{
		free(p->buf);
		p->buf = NULL;
		p->n = 0;
		p->cap = 0;
	}
}
void pa_free(pa *p, int bFreeItem)
{
	if (p == NULL)return;
	pa_clear(p, bFreeItem);
	free(p);
}

void pa_reallocate(pa* p, int cap)
{
	void** pNew;
	pNew = malloc(sizeof(void*)*cap);
	memset(pNew, 0, sizeof(void*)*cap);
	if (p->n>cap)p->n = cap;
	memcpy(pNew, p->buf, p->n * sizeof(void*));
	p->cap = cap;
	if (p->buf != NULL)free(p->buf);
	p->buf = pNew;
}
void* pa_get(pa *p, int idx/*base 0*/)
{
	if (idx<0 || idx >= p->n)return NULL;
	else return p->buf[idx];
}

void pa_set(pa *p, int idx/*base 0*/, void* item)
{
	if (idx >= p->cap)pa_reallocate(p, idx*0.2>100 ? ((int)(idx*1.2)) : idx + 100);
	p->buf[idx] = item;
	if (p->n <= idx)p->n = idx + 1;
}
//
void pa_add(pa *p, void* item)
{
	pa_set(p, p->n, item);
}

static char msg[1024 * 10];
static int gnWarning = 0;
void app_err(int nErrorCode, char *fmt, ...)
{
	va_list arg;
	/* format the message */
	va_start(arg, fmt);
	vsprintf(msg, fmt, arg);
	va_end(arg);
	printf("\n%s", msg);
#ifdef _DEBUG
	getchar();
	*((int*)NULL) = 0;
#endif
	getchar();
	exit(nErrorCode);
}

void app_warning(char *fmt, ...)
{
	va_list arg;
	gnWarning++;
	/* format the message */
	va_start(arg, fmt);
	vsprintf(msg, fmt, arg);
	va_end(arg);
	printf("\n%s", msg);
}

void app_info(char *fmt, ...)
{
	va_list arg;
	/* format the message */
	va_start(arg, fmt);
	vsprintf(msg, fmt, arg);
	va_end(arg);
	printf("\n%s", msg);
}
void app_clear_sys(app* me)
{
	int i;
	unit* pUnit;
	wind_powerPDF *pWindPowerPDF;
	wind_zone *pWindzone;

	for (i = 0; i<me->aUnits.n; i++)
	{
		pUnit = (unit*)(me->aUnits.buf[i]);
		pa_clear(&pUnit->para, 1);
		pa_clear(&pUnit->windpara, 1);
		pa_clear(&pUnit->appState, 0);
		pa_clear(&pUnit->appStateDate, 0);
		pa_clear(&pUnit->appStateDate2, 0);
		pa_clear(&pUnit->mt1, 0);
		pa_clear(&pUnit->mt2, 0);
		pa_clear(&pUnit->power, 1);
		pa_clear(&pUnit->fpower, 1);
	}
	pa_clear(&me->aUnits, 1);

	for (i = 0; i<me->aWindZone.n; i++)
	{
		pWindzone = (wind_zone*)(me->aWindZone.buf[i]);
		pa_clear(&pWindzone->para, 1);
	}
	pa_clear(&me->aWindZone, 1);

	if (&me->pmWindCorr != NULL)
	{
		mtx_free(me->pmWindCorr);
		me->pmWindCorr = NULL;
	}

	for (i = 0; i<me->aWindPowerPDF.n; i++)
	{
		pWindPowerPDF = (wind_powerPDF*)(me->aWindPowerPDF.buf[i]);
		mtx_free(pWindPowerPDF->PowerPDF);
	}
	pa_clear(&me->aWindPowerPDF, 1);
	pa_clear(&me->aSys_load, 1);
	pa_clear(&me->aRegulationPro, 1);
}

//以下一系列函数是为了得到指定日期的属性
//这些属性包括系统备用reserve，节点负荷node_para，断面属性section_para，线路属性line_para以及机组属性unit_para、指定出力unit_power以及电价属性unit_price
#define WPS_CUR_FUNCTION \
	int j; \
	for(j=0;j<WPS_CUR_PA.n;j++) \
	{ \
	if(((WPS_CUR_TYPE*)WPS_CUR_PA.buf[j])->nDate2>nDate \
	&&((WPS_CUR_TYPE*)WPS_CUR_PA.buf[j])->nDate<=nDate) \
				return (WPS_CUR_TYPE*)WPS_CUR_PA.buf[j]; \
} \
	return NULL; 

//during   

unit_para* app_cur_unit_para(unit* pUnit, int nDate)
{
#define WPS_CUR_TYPE unit_para
#define WPS_CUR_PA pUnit->para
	WPS_CUR_FUNCTION
#undef WPS_CUR_TYPE
#undef WPS_CUR_PA
}

unit_power* app_cur_unit_power(unit* pUnit, int nDate)
{
#define WPS_CUR_TYPE unit_power
#define WPS_CUR_PA pUnit->power
	WPS_CUR_FUNCTION
#undef WPS_CUR_TYPE
#undef WPS_CUR_PA
}

unit_power* app_cur_unit_fpower(unit* pUnit, int nDate)
{
#define WPS_CUR_TYPE unit_power
#define WPS_CUR_PA pUnit->fpower
	WPS_CUR_FUNCTION
#undef WPS_CUR_TYPE
#undef WPS_CUR_PA
}

wind_para* app_cur_wind_para(unit* pUnit, int nDate)
{
#define WPS_CUR_TYPE wind_para
#define WPS_CUR_PA pUnit->windpara
	WPS_CUR_FUNCTION
#undef WPS_CUR_TYPE
#undef WPS_CUR_PA
}

wind_zone_para* app_cur_windzone_para(wind_zone* pWindZone, int nDate)
{
#define WPS_CUR_TYPE wind_zone_para
#define WPS_CUR_PA pWindZone->para
	WPS_CUR_FUNCTION
#undef WPS_CUR_TYPE
#undef WPS_CUR_PA
}

wind_powerPDF* app_cur_wind_powerPDF(app* me, int nDate)
{
#define WPS_CUR_TYPE wind_powerPDF
#define WPS_CUR_PA me->aWindPowerPDF
	WPS_CUR_FUNCTION
#undef WPS_CUR_TYPE
#undef WPS_CUR_PA
}

/*获取当前系统当日负荷曲线*/
load* app_cur_load(app* me, int nDate)
{
#define WPS_CUR_TYPE load
#define WPS_CUR_PA me->aSys_load
	WPS_CUR_FUNCTION
#undef WPS_CUR_TYPE
#undef WPS_CUR_PA
}

/*获取当前日风电出力的正调峰反调峰和平出力的统计*/
regulation* app_cur_regulationpro(app* me, int nDate)
{
#define WPS_CUR_TYPE regulation
#define WPS_CUR_PA me->aRegulationPro
	WPS_CUR_FUNCTION
#undef WPS_CUR_TYPE
#undef WPS_CUR_PA
}

int app_read_para(app *me)
{
	/*ProjectID,ParaType,ItemID,BeginDay,EndDay*/
#define PT_DEBUG -1
#define PT_DEBUG_COLS 4
#define PT_DEBUG_INFO "系统调试信息"

#define PT_LOAD 2
#define PT_LOAD_COLS (POINTS_PER_DAY+5)
#define PT_LOAD_INFO "系统负荷"

#define PT_SYS 0
#define PT_SYS_COLS 18 
#define PT_SYS_INFO "系统运行参数"

#define PT_UNIT_PARA 10
#define PT_UNIT_PARA_COLS 17
#define PT_UNIT_PARA_INFO "机组参数"

#define PT_UNIT_APPSTATE 11
#define PT_UNIT_APPSTATE_COLS 6
#define PT_UNIT_APPSTATE_INFO "机组指定状态"

#define PT_UNIT_MT 13
#define PT_UNIT_MT_COLS 6
#define PT_UNIT_MT_INFO "机组检修"

#define PT_UNIT_FP 14
#define PT_UNIT_FP_COLS (POINTS_PER_DAY+5)
#define PT_UNIT_FP_INFO "机组指定出力"

#define PT_UNIT_WIND_PARA 19 
#define PT_UNIT_WIND_PARA_COLS 12
#define PT_UNIT_WIND_PARA_INFO "风电场相关信息"

#define PT_WINDZONE 50 
#define PT_WINDZONE_COLS 8
#define PT_WINDZONE_INFO "风区基本信息"

#define PT_WINDZONE_DAYAVE 51 
#define PT_WINDZONE_DAYAVE_COLS 29
#define PT_WINDZONE_DAYAVE_INFO "风区日平均风速"

#define PT_WINDZONE_MONTHAVE 52 
#define PT_WINDZONE_MONTHAVE_COLS 17
#define PT_WINDZONE_MONTHAVE_INFO "风区月平均风速"

#define PT_WINDZONE_CORR 53 
#define PT_WINDZONE_CORR_COLS 7
#define PT_WINDZONE_CORR_INFO "风区之间相关系数"


	FILE *f = NULL;
	double var[COLS];
	int nScan = 1;
	int nIdx;
	int nType, nUnitID = -1, nNodeID = -1, nLineID = -1, nSectionID = -1, nWindZoneID = -1, nID;
	int i, nDate;
	unit* pUnit;
	unit_para* pUnitPara;
	wind_para* pWindPara;
	wind_zone* pWindZone;
	wind_zone_para*pWindZonePara;
	unit_power *pUnitPower;
	load * pLoad;
	mtx* pmtemp;

	char  f_caFileFullName[255];

	int nLineCounter = 0;
	regulation *pRegulation;

	gnWarning = 0;

	pmtemp = mtx_create(0, 0, 0);

	sprintf(f_caFileFullName, "%swps_in.csv", me->caInPathName);
	printf("【提示信息】输入数据文件调用路径：%s\n", f_caFileFullName);
	f = fopen(f_caFileFullName, "r");
	if (f == NULL)
	{
		app_err(WPS_OPEN_INPUT_FILE_FAILURE,
			"无法打开输入数据文件：'WPS_in.csv'");
		return 0;
	}
	while (1)
	{
		nLineCounter++;
		if (nLineCounter == 3155)
		{
			//printf("\n测试");
		}
		nScan = fscanf(f,
			"%lf,%lf,%lf,%lf,%lf,%lf,"
			"%lf,%lf,%lf,%lf,%lf,%lf,"
			"%lf,%lf,%lf,%lf,%lf,%lf,"
			"%lf,%lf,%lf,%lf,%lf,%lf,"
			"%lf,%lf,%lf,%lf,%lf"
			, var + 0, var + 1, var + 2, var + 3, var + 4, var + 5
			, var + 6, var + 7, var + 8, var + 9, var + 10, var + 11
			, var + 12, var + 13, var + 14, var + 15, var + 16, var + 17
			, var + 18, var + 19, var + 20, var + 21, var + 22, var + 23
			, var + 24, var + 25, var + 26, var + 27, var + 28
		);
		if (nScan == EOF)break;
		if (nScan == 0)break;
		else if (nScan<2)
		{
			app_info("输入数据文件的%d行不合法可读取参数不足两个"
				, nLineCounter);
			continue;
		}
		else
		{
			/*确定数据类型*/
			nType = util_double2int(var[1]);/*消除因为浮点类型数据存储损失导致程序错误*/

											/*根据数据类型,检查数据完整性*/
			switch (nType)
			{
			case PT_SYS:
				if (nScan<PT_SYS_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_SYS_INFO, PT_SYS_COLS, nScan);
					continue;
				}
				break;
			case PT_LOAD:
				if (nScan<PT_LOAD_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_LOAD_INFO, PT_LOAD_COLS, nScan);
					continue;
				}
				break;
			case PT_UNIT_PARA:
				if (nScan<PT_UNIT_PARA_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_UNIT_PARA_INFO, PT_UNIT_PARA_COLS, nScan);
					continue;
				}
				break;
			case PT_UNIT_APPSTATE:
				if (nScan<PT_UNIT_APPSTATE_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_UNIT_APPSTATE_INFO, PT_UNIT_APPSTATE_COLS, nScan);
					continue;
				}
				break;
			case PT_UNIT_FP:
				if (nScan<PT_UNIT_FP_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_UNIT_FP_INFO, PT_UNIT_FP_COLS, nScan);
					continue;
				}
				break;
			case PT_UNIT_MT:
				if (nScan<PT_UNIT_MT_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_UNIT_MT_INFO, PT_UNIT_MT_COLS, nScan);
					continue;
				}
				break;
			case PT_UNIT_WIND_PARA:
				if (nScan<PT_UNIT_WIND_PARA_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_UNIT_WIND_PARA_INFO, PT_UNIT_WIND_PARA_COLS, nScan);
					continue;
				}
				break;
			case PT_WINDZONE:
				if (nScan<PT_WINDZONE_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_WINDZONE_INFO, PT_WINDZONE_COLS, nScan);
					continue;
				}
				break;
			case PT_WINDZONE_DAYAVE:
				if (nScan<PT_WINDZONE_DAYAVE_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_WINDZONE_DAYAVE_INFO, PT_WINDZONE_DAYAVE_COLS, nScan);
					continue;
				}
				break;
			case PT_WINDZONE_MONTHAVE:
				if (nScan<PT_WINDZONE_MONTHAVE_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_WINDZONE_MONTHAVE_INFO, PT_WINDZONE_MONTHAVE_COLS, nScan);
					continue;
				}
				break;
			case PT_WINDZONE_CORR:
				if (nScan<PT_WINDZONE_CORR_COLS)
				{
					app_warning("输入数据文件的%d行数据类型为%s，需要参数个数为%d，实际可读取参数为%d个"
						, nLineCounter, PT_WINDZONE_CORR_INFO, PT_WINDZONE_CORR_COLS, nScan);
					continue;
				}
				break;
			default:
				break;
			}

			/*要求参数文件必须是按照ID来排序的*/
			nID = util_double2int(var[2]);
			switch (nType)
			{
			case PT_UNIT_PARA:
			case PT_UNIT_APPSTATE:
			case PT_UNIT_FP:
			case PT_UNIT_MT:
			case PT_UNIT_WIND_PARA:
				if (nID<nUnitID)
				{
					app_warning("输入数据文件的%d行不合法,ID号不是由低到高排序"
						, nLineCounter);
					continue;
				}
				else if (nID == nUnitID)
				{
					pUnit = (unit*)me->aUnits.buf[me->aUnits.n - 1];
				}
				else
				{
					pUnit = malloc(sizeof(unit));
					memset(pUnit, 0, sizeof(unit));
					nUnitID = nID;
					pUnit->nId = nID;
					pa_add(&me->aUnits, pUnit);
				}
				break;
			case PT_WINDZONE:
			case PT_WINDZONE_DAYAVE:
			case PT_WINDZONE_MONTHAVE:
			case PT_WINDZONE_CORR:
				if (nID<nWindZoneID)
				{
					app_warning("输入数据文件的%d行不合法,ID号不是由低到高排序"
						, nLineCounter);
					continue;
				}
				else if (nID == nWindZoneID)
				{
					pWindZone = (wind_zone*)me->aWindZone.buf[me->aWindZone.n - 1];
				}
				else
				{
					pWindZone = malloc(sizeof(wind_zone));
					memset(pWindZone, 0, sizeof(wind_zone));
					nWindZoneID = nID;
					pWindZone->nId = nID;
					pa_add(&me->aWindZone, pWindZone);

					if (me->pmWindCorr == NULL)
					{
						me->pmWindCorr = mtx_create(1, 1, 0);
					}
					else
					{
						mtx_copy(pmtemp, me->pmWindCorr);
						mtx_resize(me->pmWindCorr, me->aWindZone.n, me->aWindZone.n);
						mtx_ini_zeros(me->pmWindCorr, me->aWindZone.n, me->aWindZone.n);
						mtx_cover(me->pmWindCorr, pmtemp);
					}
				}
				break;
			default:
				break;
			}
			/*数据检查OK了,该干活了*/
			switch (nType)
			{
			case PT_SYS:
				me->nPJID = util_double2int(var[0]);
				nIdx = 3;
				me->nStartDate = util_double2int(var[nIdx++]);
				me->nEndDate = util_double2int(var[nIdx++]);
				me->bWindFcst = util_double2int(var[nIdx++]);
				me->bCorrlation = util_double2int(var[nIdx++]);
				me->bReliability = util_double2int(var[nIdx++]);
				me->bWakeEffect = util_double2int(var[nIdx++]);
				me->bDayRhythm = util_double2int(var[nIdx++]);
				me->bSeasonRhythm = util_double2int(var[nIdx++]);
				me->nWindPowerCurve = util_double2int(var[nIdx++]);
				me->dSpeedBias = var[nIdx++];
				me->nErrorType = util_double2int(var[nIdx++]);
				me->nErrorCoherence = util_double2int(var[nIdx++]);
				me->nRandomType = util_double2int(var[nIdx++]);
				me->nRandomSeed = util_double2int(var[nIdx++]);
				me->nSimulateTimes = util_double2int(var[nIdx++]);
				me->bAppPower = util_double2int(var[nIdx++]);
				me->dPeakValleyRegulationLimit = WPS_POS_NEG_REGULATION;
				break;
			case PT_LOAD:
				pLoad = (load*)malloc(sizeof(load));
				memset(pLoad, 0, sizeof(load));
				pa_add(&me->aSys_load, pLoad);
				nIdx = 3;
				pLoad->nDate = util_double2int(var[nIdx++]);				//起始时间
				pLoad->nDate2 = util_double2int(var[nIdx++]);				//终止时间
				memcpy(pLoad->dLoad, var + nIdx, sizeof(double)*POINTS_PER_DAY);
				break;
			case PT_UNIT_APPSTATE:
				nIdx = 3;
				pa_add(&pUnit->appStateDate, (void*)(util_double2int(var[nIdx++])));
				pa_add(&pUnit->appStateDate2, (void*)(util_double2int(var[nIdx++])));
				pa_add(&pUnit->appState, (void*)(util_double2int(var[nIdx++])));
				break;
			case PT_UNIT_FP:
				pUnitPower = (unit_power*)malloc(sizeof(unit_power));
				memset(pUnitPower, 0, sizeof(unit_power));
				pa_add(&pUnit->power, pUnitPower);
				nIdx = 3;
				pUnitPower->nDate = util_double2int(var[nIdx++]);
				pUnitPower->nDate2 = util_double2int(var[nIdx++]);
				pUnitPower->nFixedPower = 1;
				memcpy(pUnitPower->dPower, var + nIdx, sizeof(double)*POINTS_PER_DAY);
				break;
			case PT_UNIT_PARA:
				pUnitPara = malloc(sizeof(unit_para));
				memset(pUnitPara, 0, sizeof(unit_para));
				pa_add(&pUnit->para, pUnitPara);
				nIdx = 3;
				pUnitPara->nDate = util_double2int(var[nIdx++]);
				pUnitPara->nDate2 = util_double2int(var[nIdx++]);
				pUnitPara->nNodeId = util_double2int(var[nIdx++]);
				pUnitPara->nPlantId = util_double2int(var[nIdx++]);
				pUnitPara->nCompanyID = util_double2int(var[nIdx++]);
				pUnitPara->nUnitTypeIn = util_double2int(var[nIdx++]);
				pUnitPara->dPowerMax = var[nIdx++];
				pUnitPara->dPowerMin = var[nIdx++];
				if (pUnitPara->dPowerMax<pUnitPara->dPowerMin - TINY)
				{
					app_warning("输入数据文件的%d行不合法,机组最小出力大于机组容量", nLineCounter);
					continue;
				}
				pUnitPara->dEnergyCapa = var[nIdx++];
				pUnitPara->dConvertEfft = var[nIdx++];
				pUnitPara->dUpDownCost = var[nIdx++];
				pUnitPara->dFOR = var[nIdx++];
				pUnitPara->dSelfUse = var[nIdx++];
				pUnitPara->dVariableCost = var[nIdx++];
				pUnitPara->dVariableCost2 = var[nIdx++];
				pUnitPara->dCoal = var[nIdx++];
				pUnitPara->dCoal2 = var[nIdx++];
				break;
			case PT_UNIT_WIND_PARA:
				pWindPara = malloc(sizeof(wind_para));
				memset(pWindPara, 0, sizeof(wind_para));
				pa_add(&pUnit->windpara, pWindPara);
				nIdx = 3;
				pWindPara->nDate = util_double2int(var[nIdx++]);
				pWindPara->nDate2 = util_double2int(var[nIdx++]);
				pWindPara->nWindZone = util_double2int(var[nIdx++]);
				pWindPara->dCutinSP = var[nIdx++];
				pWindPara->dRatedSP = var[nIdx++];
				pWindPara->dCutoutSP = var[nIdx++];
				pWindPara->dWakeeffect = var[nIdx++];
				pWindPara->dAvailableRate = var[nIdx++];
				pWindPara->dForecastError = var[nIdx++];
				break;
			case PT_UNIT_MT:
				nIdx = 3;
				pa_add(&pUnit->mt1, (void*)(util_double2int(var[nIdx++])));
				//检修计划里面日期系统中，检修日包括检修截止的那天，与运行模拟里面的日期系统不同，因此这里截止日期要加一天，zhangning
				pa_add(&pUnit->mt2, (void*)(idate_next_day(util_double2int(var[nIdx++]))));
				break;
			case PT_WINDZONE:
				pWindZonePara = malloc(sizeof(wind_zone_para));
				memset(pWindZonePara, 0, sizeof(wind_zone_para));
				pa_add(&pWindZone->para, pWindZonePara);
				nIdx = 3;
				pWindZonePara->nDate = util_double2int(var[nIdx++]);
				pWindZonePara->nDate2 = util_double2int(var[nIdx++]);
				pWindZonePara->dWeibull_c = var[nIdx++];
				pWindZonePara->dWeibull_k = var[nIdx++];
				pWindZonePara->dtheta = var[nIdx++];
				break;
			case PT_WINDZONE_DAYAVE:
				nIdx = 3;
				nDate = util_double2int(var[nIdx++]);
				pWindZonePara = app_cur_windzone_para(pWindZone, nDate);
				nIdx = 5;
				memcpy(pWindZonePara->nDayAve, var + nIdx, sizeof(double)*POINTS_PER_DAY);
				break;
			case PT_WINDZONE_MONTHAVE:
				nIdx = 3;
				nDate = util_double2int(var[nIdx++]);
				pWindZonePara = app_cur_windzone_para(pWindZone, nDate);
				nIdx = 5;
				memcpy(pWindZonePara->nMonthAve, var + nIdx, sizeof(double) * 12);
				break;
			case PT_WINDZONE_CORR:
				nIdx = 5;
				for (i = 1; i <= me->aWindZone.n; i++)
				{
					pWindZone = (wind_zone*)me->aWindZone.buf[i - 1];
					if (pWindZone->nId == util_double2int(var[nIdx]))
					{
						mtx_set_element(me->pmWindCorr, me->aWindZone.n, i, var[nIdx + 1]);
						mtx_set_element(me->pmWindCorr, i, me->aWindZone.n, var[nIdx + 1]);
					}
				}
				break;
			default:
				break;
			}
		}
	}
	for (i = idate_first_day_of_year(me->nStartDate); i<me->nEndDate; i = idate_next_year(i))
	{
		pRegulation = malloc(sizeof(regulation));
		memset(pRegulation, 0, sizeof(regulation));
		pRegulation->nDate = i;
		pRegulation->nDate2 = idate_next_year(i);
		pa_add(&me->aRegulationPro, pRegulation);
	}

	mtx_free(pmtemp);
	return 1;
}


void   app_read_calDir(app *me)
{
	FILE *fp;
	int  f_nNum = 0;
	char f_TxtDir[255];

	/*获取当前exe计算目录，实际应用中可能是在另外的目录中（比如D盘）来调用这个EXE的，
	而RECC_cal_dir.txt是跟exe放在同一个目录下的，这种情况下如果直接读就会默认读D盘，导致程序错误。*/
	strcpy(me->sCursDir, me->caExePathName);
	if (strrchr(me->sCursDir, '\\') == NULL)//strrchr()函数查找字符在制定字符串中从后面开始的第一次的出现位置
	{
		sprintf(f_TxtDir, "WPS_cal_dir.txt");
		printf("【提示信息】数据文件目录：\n");
	}
	else
	{
		strrchr(me->sCursDir, '\\')[0] = 0;
		printf("【提示信息】数据文件目录：%s\n", me->sCursDir);
		sprintf(f_TxtDir, "%s\\WPS_cal_dir.txt", me->sCursDir);
	}
	//printf("读取文件%s\n",f_CursDir);

	fp = fopen(f_TxtDir, "r");
	while (!feof(fp) && f_nNum<WPS_MAX_PROJECT)
	{
		if (fscanf(fp, "%s", me->sDir[f_nNum]) >= 1)
		{
			f_nNum++;
		}
	}
	me->nDir = f_nNum;
	fclose(fp);

}


int wpsmain(int argc, char*argv[])
{
	app sys;
	int j;
	char   f_caFullFileName[255];
#ifdef WIN32
	/*如果是Windows系统，将程序优先级设置为最低*/
	SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS);
#endif
	sys.t_start = clock();
	memset(&sys, 0, sizeof(app));
	strcpy(sys.caExePathName, argv[0]);
	printf("【提示信息】程序调用路径：%s\n", sys.caExePathName);
	/*读入计算目录*/
	app_read_calDir(&sys);

	for (j = 0; j<sys.nDir; j++)
	{
		sprintf(sys.caInPathName, "%s\\%s", sys.sCursDir, sys.sDir[j]);
		sprintf(sys.caOutPathName, "%s\\%s", sys.sCursDir, sys.sDir[j]);
		mkdir(sys.caOutPathName);

		/*确定黑窗口信息输出目录*/
		if (strrchr(sys.caExePathName, '\\') == NULL)//strrchr()函数查找字符在制定字符串中从后面开始的第一次的出现位置
		{
			printf("【提示信息】当前数据文件路径：\\%s\n", sys.sDir[j]);
		}
		else
		{
			printf("【提示信息】当前数据文件路径：%s\\%s\n", sys.sCursDir, sys.sDir[j]);
		}

		/*读入基础数据*/
		app_read_para(&sys);

		if (gnWarning>0)
		{
			app_err(WPS_EXIT_PARAERROR, "输入参数有误");
		}

		sprintf(f_caFullFileName, "%swps_out.csv", sys.caOutPathName);
		sys.fRet = fopen(f_caFullFileName, "w+");
		if (sys.fRet == NULL)
		{
			app_err(WPS_OPEN_OUTPUT_FILE_FAILURE, "无法写入数据文件：'wps_out.csv'");
			return;
		}
		for (sys.nTimes = 1; sys.nTimes <= sys.nSimulateTimes; sys.nTimes++)
		{
			printf("[提示信息]正在进行第%d场景(共%d场景)风电场出力模拟中……\n", sys.nTimes, sys.nSimulateTimes);

			/*生成风电出力*/
			app_generate_wind_power(&sys);

			/*输出模拟结果*/
			app_out_wind_power(&sys);

			app_out_wind_power_stat(&sys);

			/*输出风电预测出力模拟结果*/
			if (sys.bWindFcst == 1)
			{
				app_out_wind_fpower(&sys);
				app_out_wind_ferror(&sys);
				app_out_wind_fpower_stat(&sys);
			}
		}

		/*释放资源*/
		app_clear_sys(&sys);
		/*内存泄露检测*/
		lib_check_mem_leak();
		sys.t_finish = clock();
	}
	printf("\n计算程序所用时间：%.10g秒", ((double)(sys.t_finish - sys.t_start)) / CLOCKS_PER_SEC);
	//getchar();
	return WPS_EXIT_SUCCESS;
}



BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

