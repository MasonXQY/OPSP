#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include <math.h>
#include <stdarg.h>
#include <LIMITS.H>
#include <time.h>
#include "globle.h"
#include "gConstant.h"
#include "func.h"
#pragma warning(disable: 4996)

static void fault(char *fmt, ...);

/*
该函数比较重要，一般情况下请勿随便修改
问题：浮点数值为2，但经过一系列运算后可能结果为1.999999999999999999999999
使用类型强制转换为int，结果为1
*/
int util_double2int(double d)
{
	if(d>=0&&d<INT_MAX)
		return (int)(d+.5);
	else if(d<=0&&d>INT_MIN)
		return (int)(d-.5);
	else
	{
		return (int)d;
	}
}
/*HOOK printf*/
static FILE* fInfo=NULL;
static void WPS_printf(char *fmt, ...)
{
      va_list arg;
      char msg[4095+1];
      /* format the message */
      va_start(arg, fmt);
      vsprintf(msg, fmt, arg);
      va_end(arg);
	  fprintf(stdout,"%s",msg);
	  if(fInfo==NULL)
	  {
		  fInfo=fopen("WPS_RunInfo.txt","w+");
	  }
	  if(fInfo!=NULL)
	  {
		  fprintf(fInfo,"%s",msg);
		  fflush(fInfo);
	  }
}
#define printf WPS_printf

static void fault(char *fmt, ...)
{
      va_list arg;
      char msg[4095+1];
      /* format the message */
      va_start(arg, fmt);
      vsprintf(msg, fmt, arg);
      va_end(arg);
      /* send the message to the standard output */
      fprintf(stdout, "%s\n", msg);
	  /* terminate program execution */
#ifdef _DEBUG
	  *((int*)NULL)=0;
#endif
      exit(EXIT_FAILURE);
      /* no return */
}

//during
void mtx_check(mtx* temp)
{
	double *data;
	int i;
	data=mtx_get_buf_read(temp);
	for (i=0;i<temp->m*temp->n;i++)
		if(*(data+i)<TINY) *(data+i)=0;
}

#define EX_LEN 5
static double ex[EX_LEN];
void app_out_ex_clear()
{
	memset(ex,0,sizeof(double)*EX_LEN);
}
void app_out_ex_set(int nID,double var)
{
	if (nID<0||nID>EX_LEN)
		app_err(WPS_EXIT_PARAERROR,"参数错误");
	ex[nID-1]=var;
}
void app_out_element(FILE* f,double var)
{
	fprintf(f,",%.10g",var);
}
void app_out_array(app* me,int nPJID,int nItemID,int nValueType,int nDate,int nNextDate,double var[],int bIgnoreAllZeros)
{
	int i=0,j=0;
	double dSumPos=0,dSumNeg=0;
	int nMaxIdx=0,nMinIdx=0;
	//int nRunLevels;
	//int nRealNextDate=idate_next_month(nDate);/*扩展为每个月1号肯定有数据*/
	int nRealNextDate=nNextDate;

	for(j=0;j<POINTS_PER_DAY;j++)
	{
		if(var[j]>0)dSumPos+=var[j];else dSumNeg+=var[j];
		if(var[j]>var[nMaxIdx])nMaxIdx=j;
		if(var[j]<var[nMinIdx])nMinIdx=j;
	}

	if(bIgnoreAllZeros!=0)if(dSumPos-dSumNeg<1e-3)return;
	if(nRealNextDate>=nNextDate)nRealNextDate=nNextDate;
	fprintf(me->fRet,"%d,%d,%d,%d,%d,%d,%d,%d",nPJID,nValueType,nItemID,nDate
		,nRealNextDate,idate_days_between(nDate,nRealNextDate),nDate/10000,(nDate/100)%100);
	for(j=0;j<POINTS_PER_DAY;j++)app_out_element(me->fRet,var[j]);
	app_out_element(me->fRet,dSumPos+dSumNeg);
	app_out_element(me->fRet,dSumPos);
	app_out_element(me->fRet,dSumNeg);
	app_out_element(me->fRet,var[nMaxIdx]);
	app_out_element(me->fRet,var[nMinIdx]);
	for(j=0;j<EX_LEN;j++)app_out_element(me->fRet,ex[j]);
	fprintf(me->fRet,"\n");

	if(nRealNextDate<nNextDate)
		app_out_array(me,nPJID,nItemID,nValueType,nRealNextDate,nNextDate,var,bIgnoreAllZeros);
}

void app_stat_wind_power( double *dapower, int nCount, double dPowerMax, double *daPDF, double *daCDF, double *dminoutput,double *dmaxoutput,double *dcapacityfactor, double dconfidence_min, double dconfidence_max)
{
	int i,j;
	double dstep=dPowerMax/WPS_DISTRIBUTION_INTERVAL;
	double *dapowertemp;
	double dtemp;
	double dmean=0;

	dapowertemp=(double*)malloc(sizeof(double)*nCount);
	memset(dapowertemp,0,sizeof(double)*nCount);
	memset(daPDF,0,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
	memset(daCDF,0,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);

	//统计PDF
	for(i=0;i<nCount;i++)	
	{
		j=util_double2int(dapower[i]/dstep);
		if (j>=WPS_DISTRIBUTION_INTERVAL)
		{
			j=WPS_DISTRIBUTION_INTERVAL-1;
		}
		daPDF[j]=daPDF[j]+1.0/nCount;
		dmean+=dapower[i]/nCount;
	}
	//统计CDF
	daCDF[0]=daPDF[0];
	for(i=1;i<WPS_DISTRIBUTION_INTERVAL;i++)	
	{
		daCDF[i]=daCDF[i-1]+daPDF[i];
	}
	//对出力进行冒泡排序，从小到大
	memcpy(dapowertemp,dapower,sizeof(double)*nCount);
	for(i=0;i<nCount-1;i++)	
	{
		for(j=0;j<nCount-1-i;j++)	
		{
			if (dapowertemp[j]>dapowertemp[j+1])
			{
				dtemp=dapowertemp[j];
				dapowertemp[j]=dapowertemp[j+1];
				dapowertemp[j+1]=dtemp;
			}
		}
	}
	*dminoutput=dapowertemp[util_double2int((1-dconfidence_min)*nCount)];
	j=util_double2int(dconfidence_max*nCount);	
	if (j==nCount) j--;
	*dmaxoutput=dapowertemp[j];
	if (*dmaxoutput>=dPowerMax)
		*dmaxoutput=dPowerMax-TINY;
	*dcapacityfactor=dmean;
	free(dapowertemp);
}

void app_out_wind_power_stat( app* me ) 
{
	int i,j,nCount,nCount2;
	int nDate,nYear,ntemp;
	int nPosDay=0,nFlatDay=0,nNegDay=0;
	unit *pUnit;
	unit_para *pUnitPara,*pUnitParatemp;
	wind_para* pWindPara,*pWindParatemp;
	unit_power* pUnitPower;
	double dapower[POINTS_PER_DAY*DAYS_PER_YEAR]={0};
	double dapowerall[POINTS_PER_DAY*DAYS_PER_YEAR]={0};
	int napowerall[POINTS_PER_DAY*DAYS_PER_YEAR]={0};
	double var[POINTS_PER_DAY]={0};
	int nBeginYear=0,nEndYear=0;
	double daPDF[WPS_DISTRIBUTION_INTERVAL]={0};
	double daCDF[WPS_DISTRIBUTION_INTERVAL]={0};
	double dminoutput,dmaxoutput,dmeanoutput,dPowerMax,dPowerMaxall,dmaxoutputSigma;
	double dPosRegulation[POINTS_PER_DAY]={0};
	double dNegRegulation[POINTS_PER_DAY]={0};
	double dFlatRegulation[POINTS_PER_DAY]={0};
	double dAverageLoad[POINTS_PER_DAY]={0};
	double dMaxHourlyOutput[POINTS_PER_DAY]={0};
	double dMinHourlyOutput[POINTS_PER_DAY]={0};
	double dAverageHourlyOutput[POINTS_PER_DAY]={0};
	double *daLoad;
	double dPosRegulationPro,dNegRegulationPro,dFlatRegulationPro;
	load* pLoad;
	regulation *pRegulationPro;
	daLoad=(double*)malloc(sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);


	//由于CDF和PDF都需要通过var输出，所以这里规定其长度不能大于var
	if (WPS_DISTRIBUTION_INTERVAL>POINTS_PER_DAY-4)
	{
		app_err(WPS_EXIT_PARAERROR,"输入参数有误");
	}

	nBeginYear=(me->nStartDate/10000);
	if (me->nEndDate==idate_first_day_of_year(me->nEndDate))
	{
		nEndYear=(me->nEndDate/10000)-1;
	}
	else
	{
		nEndYear=(me->nEndDate/10000);
	}
	for (nYear=nBeginYear;nYear<=nEndYear;nYear++)
	{
		dPowerMaxall=0;dmaxoutputSigma=0;
		memset(dapowerall,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
		memset(napowerall,0,sizeof(int)*POINTS_PER_DAY*DAYS_PER_YEAR);

		//输出风电场典型出力
		memset(dPosRegulation,0,sizeof(double)*POINTS_PER_DAY);
		memset(dNegRegulation,0,sizeof(double)*POINTS_PER_DAY);
		memset(dFlatRegulation,0,sizeof(double)*POINTS_PER_DAY);
		memset(dAverageLoad,0,sizeof(double)*POINTS_PER_DAY);
		memset(daLoad,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);

		for(i=0;i<me->aUnits.n;i++)
		{
			memset(dapower,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
			pUnit=(unit*)(me->aUnits.buf[i]);
			nCount=0;dPowerMax=0;
			for (nDate=nYear*10000+101;nDate<nYear*10000+10101;nDate=idate_next_day(nDate))
			{
				if (nDate<me->nStartDate||nDate>=me->nEndDate) continue;
				pUnitPara=app_cur_unit_para(pUnit,nDate);
				if(pUnitPara==NULL) continue;
				if(pUnitPara->nUnitTypeIn!=11) continue;		
				pWindPara=app_cur_wind_para(pUnit,nDate);
				if (pWindPara==NULL) continue;	
				if (dPowerMax<pUnitPara->dPowerMax)
				{
					dPowerMax=pUnitPara->dPowerMax;
				}
				pUnitParatemp=pUnitPara;pWindParatemp=pWindPara;
				pUnitPower=app_cur_unit_power(pUnit,nDate);
				if (pUnitPower==NULL) continue;	
				ntemp=idate_day_idx_of_year(nDate);
				for(j=0;j<POINTS_PER_DAY;j++)	
				{
					dapower[nCount]=pUnitPower->dPower[j];
					nCount++;
					dapowerall[(ntemp-1)*POINTS_PER_DAY+j]+=pUnitPower->dPower[j];
					napowerall[(ntemp-1)*POINTS_PER_DAY+j]=1;
				}
			}
			//开始统计风电特性
			if (nCount>0)
			{
				dPowerMaxall+=dPowerMax;
				app_stat_wind_power(dapower,nCount, dPowerMax,daPDF,daCDF,&dminoutput,&dmaxoutput,&dmeanoutput,WPS_MIN_OUTPUT_CONFIDENCE,WPS_MAX_OUTPUT_CONFIDENCE);
				dmaxoutputSigma+=dmaxoutput;
				app_out_ex_clear();
				app_out_ex_set(1,me->nTimes);
				app_out_ex_set(2,dPowerMax);
				app_out_ex_set(3,pUnitParatemp->nNodeId);
				app_out_ex_set(4,pWindParatemp->nWindZone);
				memset(var,0,sizeof(double)*POINTS_PER_DAY);
				var[POINTS_PER_DAY-3]=dminoutput;
				var[POINTS_PER_DAY-2]=dmaxoutput;
				var[POINTS_PER_DAY-1]=dmeanoutput;
				memcpy(var,daPDF,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
				app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_STAT_PDF,nYear*10000+101,nYear*10000+10101,var,0);
				memcpy(var,daCDF,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
				app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_STAT_CDF,nYear*10000+101,nYear*10000+10101,var,0);
			}
		}
		nCount=0;
		memset(dapower,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
		for(j=0;j<POINTS_PER_DAY*DAYS_PER_YEAR;j++)	
		{
			if (napowerall[j]==1)
			{
				dapower[nCount]=dapowerall[j];
				nCount++;
			}
		}
		if (nCount>0)
		{
			app_stat_wind_power(dapower,nCount, dPowerMaxall,daPDF,daCDF,&dminoutput,&dmaxoutput,&dmeanoutput,WPS_MIN_OUTPUT_CONFIDENCE,WPS_MAX_OUTPUT_CONFIDENCE);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			var[POINTS_PER_DAY-4]=dmaxoutput/dmaxoutputSigma;
			var[POINTS_PER_DAY-3]=dminoutput;
			var[POINTS_PER_DAY-2]=dmaxoutput;
			var[POINTS_PER_DAY-1]=dmeanoutput;
			memcpy(var,daPDF,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
			app_out_array(me,me->nPJID,0,RT_WIND_ALL_STAT_PDF,nYear*10000+101,nYear*10000+10101,var,0);
			memcpy(var,daCDF,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
			app_out_array(me,me->nPJID,0,RT_WIND_ALL_STAT_CDF,nYear*10000+101,nYear*10000+10101,var,0);
			//输出风电场典型出力
			memset(dPosRegulation,0,sizeof(double)*POINTS_PER_DAY);
			memset(dNegRegulation,0,sizeof(double)*POINTS_PER_DAY);
			memset(dFlatRegulation,0,sizeof(double)*POINTS_PER_DAY);
			memset(dAverageLoad,0,sizeof(double)*POINTS_PER_DAY);
			memset(dMaxHourlyOutput,0,sizeof(double)*POINTS_PER_DAY);
			memset(dMinHourlyOutput,0,sizeof(double)*POINTS_PER_DAY);
			memset(dAverageHourlyOutput,0,sizeof(double)*POINTS_PER_DAY);
			memset(daLoad,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
			nCount2=0;
			for (nDate=nYear*10000+101;nDate<nYear*10000+10101;nDate=idate_next_day(nDate))
			{
				pLoad=app_cur_load(me,nDate);
				for(j=0;j<POINTS_PER_DAY;j++)	
				{
					daLoad[nCount2+j]=pLoad->dLoad[j];
					dAverageLoad[j]+=pLoad->dLoad[j]/nCount*POINTS_PER_DAY;
				}
				for(i=0;i<me->aUnits.n;i++)
				{
					pUnit=(unit*)(me->aUnits.buf[i]);
					pUnitPara=app_cur_unit_para(pUnit,nDate);
					if (pUnitPara!=NULL&&pUnitPara->nUnitTypeIn==11) continue;	
					pUnitPower=app_cur_unit_power(pUnit,nDate);
					if (pUnitPower==NULL) continue;	
					for(j=0;j<POINTS_PER_DAY;j++)	
					{
						if (pUnitPower->dPower[j]<0)
						{
							daLoad[nCount2+j]-=pUnitPower->dPower[j];
							dAverageLoad[j]-=pUnitPower->dPower[j]/nCount*POINTS_PER_DAY;
						}
					}
				}
				nCount2+=POINTS_PER_DAY;
			}

			dPosRegulationPro=0;dNegRegulationPro=0;dFlatRegulationPro=0;
			pRegulationPro=app_cur_regulationpro(me,nYear*10000+101);
			app_stat_wind_typical_curve(me, daLoad,dapower,nCount,dPosRegulation,dNegRegulation,dFlatRegulation,&dPosRegulationPro,&dNegRegulationPro,&dFlatRegulationPro);
			dPosRegulationPro=0;
			dNegRegulationPro=0;
			dFlatRegulationPro=0;
			for(i=0;i<12;i++)
			{
				nPosDay+=pRegulationPro->nPosRegulationdays[i];
				nFlatDay+=pRegulationPro->nFlatRegulationdays[i];
				nNegDay+=pRegulationPro->nNegRegulationdays[i];
				dPosRegulationPro+=pRegulationPro->nPosRegulationdays[i]*pRegulationPro->dPosRegulationPro[i];
				dFlatRegulationPro+=pRegulationPro->nFlatRegulationdays[i]*pRegulationPro->dFlatRegulationPro[i];
				dNegRegulationPro+=pRegulationPro->nNegRegulationdays[i]*pRegulationPro->dNegRegulationPro[i];
			}
			dPosRegulationPro/=(nPosDay+TINY);
			dFlatRegulationPro/=(nFlatDay+TINY);
			dNegRegulationPro/=(nNegDay+TINY);
			//输出模拟典型出力
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,dPosRegulationPro);
			app_out_ex_set(4,0);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			memcpy(var,dPosRegulation,sizeof(double)*POINTS_PER_DAY);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_SIMULATE,nYear*10000+101,nYear*10000+10101,var,0);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,dFlatRegulationPro);
			app_out_ex_set(4,1);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			memcpy(var,dFlatRegulation,sizeof(double)*POINTS_PER_DAY);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_SIMULATE,nYear*10000+101,nYear*10000+10101,var,0);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,dNegRegulationPro);
			app_out_ex_set(4,2);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			memcpy(var,dNegRegulation,sizeof(double)*POINTS_PER_DAY);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_SIMULATE,nYear*10000+101,nYear*10000+10101,var,0);
			//输出理论典型出力
			app_stat_wind_daily_curve(me, dAverageLoad,dapower,nCount,dPowerMaxall,dPosRegulation,dNegRegulation,dMaxHourlyOutput,dMinHourlyOutput,dAverageHourlyOutput);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,dPosRegulationPro);
			app_out_ex_set(4,0);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			memcpy(var,dPosRegulation,sizeof(double)*POINTS_PER_DAY);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_THEORETICAL,nYear*10000+101,nYear*10000+10101,var,0);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,dFlatRegulationPro);
			app_out_ex_set(4,1);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			memcpy(var,dAverageHourlyOutput,sizeof(double)*POINTS_PER_DAY);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_THEORETICAL,nYear*10000+101,nYear*10000+10101,var,0);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,dNegRegulationPro);
			app_out_ex_set(4,2);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			memcpy(var,dNegRegulation,sizeof(double)*POINTS_PER_DAY);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_THEORETICAL,nYear*10000+101,nYear*10000+10101,var,0);
		}
	}
	free(daLoad);
}

//该函数与app_out_wind_power_stat相同，只是按月统计输出，需要在mtx结构体里读取，模拟无效的mtx也计入统计中。
void app_out_wind_power_stat_month( app* me, mtx* pmWindPower[WPS_SIMULATION_TIMES],int nStartDate, int nEndDate) 
{
	int i,j,nCount=0,k,l,nWindFarm,nCount2;
	int nDate,nmonth;
	unit *pUnit;
	unit_para *pUnitPara;
	wind_para* pWindPara;
	double *dapower;
	double *dapowerall;
	int *napowerall;
	double var[POINTS_PER_DAY]={0};
	int nBeginYear=0,nEndYear=0;
	double daPDF[WPS_DISTRIBUTION_INTERVAL]={0};
	double daCDF[WPS_DISTRIBUTION_INTERVAL]={0};
	double dminoutput,dmaxoutput,dmeanoutput,dPowerMax,dPowerMaxall,dmaxoutputSigma;
	double dPosRegulation[POINTS_PER_DAY]={0};
	double dNegRegulation[POINTS_PER_DAY]={0};
	double dFlatRegulation[POINTS_PER_DAY]={0};
	double dAverageLoad[POINTS_PER_DAY]={0};
	double dMaxHourlyOutput[POINTS_PER_DAY]={0};
	double dMinHourlyOutput[POINTS_PER_DAY]={0};
	double dAverageHourlyOutput[POINTS_PER_DAY]={0};
	double *daLoad;
	double dPosRegulationPro,dNegRegulationPro,dFlatRegulationPro;
	load* pLoad;
	unit_power* pUnitPower;
	regulation * pRegulationPro;

	daLoad=(double*)malloc(sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	dapower=(double*)malloc(sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	dapowerall=(double*)malloc(sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	napowerall=(int*)malloc(sizeof(int)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);

	//由于CDF和PDF都需要通过var输出，所以这里规定其长度不能大于var
	if (WPS_DISTRIBUTION_INTERVAL>POINTS_PER_DAY-4)
	{
		app_err(WPS_EXIT_PARAERROR,"输入参数有误");
	}

	nmonth= idate_month_idx_of_year(nStartDate);
	nDate=nStartDate;

	dPowerMaxall=0;dmaxoutputSigma=0;
	memset(dapowerall,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	memset(napowerall,0,sizeof(int)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	nWindFarm=0;
	for(i=0;i<me->aUnits.n;i++)
	{
		memset(dapower,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
		pUnit=(unit*)(me->aUnits.buf[i]);
		pUnitPara=app_cur_unit_para(pUnit,nDate);
		if (pUnitPara==NULL) continue;
		if(pUnitPara->nUnitTypeIn!=11)continue;		
		pWindPara=app_cur_wind_para(pUnit,nDate);
		if (pWindPara==NULL) continue;	
		nCount=0;
		dPowerMax=pUnitPara->dPowerMax;
		nWindFarm++;
		for (l=0;l<WPS_SIMULATION_TIMES;l++)
		{
			for (k=1;k<= pmWindPower[l]->n;k++)
			{
				dapower[nCount]=mtx_get_element(pmWindPower[l],nWindFarm,k);
				dapowerall[nCount]+=mtx_get_element(pmWindPower[l],nWindFarm,k);
				napowerall[nCount]=1;
				nCount++;
			}
		}
		//开始统计风电特性
		if (nCount>0)
		{
			dPowerMaxall+=dPowerMax;
			app_stat_wind_power(dapower,nCount, dPowerMax,daPDF,daCDF,&dminoutput,&dmaxoutput,&dmeanoutput,WPS_MIN_OUTPUT_CONFIDENCE,WPS_MAX_OUTPUT_CONFIDENCE);
			dmaxoutputSigma+=dmaxoutput;
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMax);
			app_out_ex_set(3,pUnitPara->nNodeId);
			app_out_ex_set(4,pWindPara->nWindZone);
			app_out_ex_set(5,nmonth);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			var[POINTS_PER_DAY-3]=dminoutput;
			var[POINTS_PER_DAY-2]=dmaxoutput;
			var[POINTS_PER_DAY-1]=dmeanoutput;
			memcpy(var,daPDF,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_STAT_PDF,nStartDate,nEndDate,var,0);
			memcpy(var,daCDF,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_STAT_CDF,nStartDate,nEndDate,var,0);
		}
	}
	
	memset(dapower,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	for(j=0;j<nCount;j++)	
	{
		if (napowerall[j]==1)
		{
			dapower[j]=dapowerall[j];
		}
	}
	if (nCount>0)
	{
		//输出风电场总体统计参数
		app_stat_wind_power(dapower,nCount, dPowerMaxall,daPDF,daCDF,&dminoutput,&dmaxoutput,&dmeanoutput,WPS_MIN_OUTPUT_CONFIDENCE,WPS_MAX_OUTPUT_CONFIDENCE);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		var[POINTS_PER_DAY-4]=dmaxoutput/dmaxoutputSigma;
		var[POINTS_PER_DAY-3]=dminoutput;
		var[POINTS_PER_DAY-2]=dmaxoutput;
		var[POINTS_PER_DAY-1]=dmeanoutput;
		memcpy(var,daPDF,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
		app_out_array(me,me->nPJID,0,RT_WIND_ALL_STAT_PDF,nStartDate,nEndDate,var,0);
		memcpy(var,daCDF,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
		app_out_array(me,me->nPJID,0,RT_WIND_ALL_STAT_CDF,nStartDate,nEndDate,var,0);
		//输出风电场典型出力
		memset(dPosRegulation,0,sizeof(double)*POINTS_PER_DAY);
		memset(dNegRegulation,0,sizeof(double)*POINTS_PER_DAY);
		memset(dFlatRegulation,0,sizeof(double)*POINTS_PER_DAY);
		memset(dAverageLoad,0,sizeof(double)*POINTS_PER_DAY);
		memset(dMaxHourlyOutput,0,sizeof(double)*POINTS_PER_DAY);
		memset(dMinHourlyOutput,0,sizeof(double)*POINTS_PER_DAY);
		memset(dAverageHourlyOutput,0,sizeof(double)*POINTS_PER_DAY);
		memset(daLoad,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
		nCount2=0;
		for (l=0;l<WPS_SIMULATION_TIMES;l++)
		{
			for (nDate=nStartDate;nDate<nEndDate;nDate=idate_next_day(nDate))
			{
				pLoad=app_cur_load(me,nDate);
				for(j=0;j<POINTS_PER_DAY;j++)	
				{
					daLoad[nCount2+j]=pLoad->dLoad[j];
					dAverageLoad[j]+=pLoad->dLoad[j]/nCount*POINTS_PER_DAY;
				}
				for(i=0;i<me->aUnits.n;i++)
				{
					pUnit=(unit*)(me->aUnits.buf[i]);
					pUnitPara=app_cur_unit_para(pUnit,nDate);
					if (pUnitPara!=NULL&&pUnitPara->nUnitTypeIn==11) continue;	
					pUnitPower=app_cur_unit_power(pUnit,nDate);
					if (pUnitPower==NULL) continue;	
					for(j=0;j<POINTS_PER_DAY;j++)	
					{
						if (pUnitPower->dPower[j]<0)
						{
							daLoad[nCount2+j]-=pUnitPower->dPower[j];
							dAverageLoad[j]-=pUnitPower->dPower[j]/nCount*POINTS_PER_DAY;
						}
					}
				}
				nCount2+=POINTS_PER_DAY;
			}
		}

		dPosRegulationPro=0;dNegRegulationPro=0;dFlatRegulationPro=0;
		pRegulationPro=app_cur_regulationpro(me,nStartDate);
		app_stat_wind_typical_curve(me, daLoad,dapower,nCount,dPosRegulation,dNegRegulation,dFlatRegulation,&dPosRegulationPro,&dNegRegulationPro,&dFlatRegulationPro);
		pRegulationPro->nPosRegulationdays[nmonth-1]=idate_days_between(nStartDate,nEndDate);
		pRegulationPro->nFlatRegulationdays[nmonth-1]=idate_days_between(nStartDate,nEndDate);
		pRegulationPro->nNegRegulationdays[nmonth-1]=idate_days_between(nStartDate,nEndDate);
		pRegulationPro->dPosRegulationPro[nmonth-1]=dPosRegulationPro;
		pRegulationPro->dFlatRegulationPro[nmonth-1]=dFlatRegulationPro;
		pRegulationPro->dNegRegulationPro[nmonth-1]=dNegRegulationPro;
		//printf("测试\n");
		//输出模拟典型出力
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,dPosRegulationPro);
		app_out_ex_set(4,0);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		memcpy(var,dPosRegulation,sizeof(double)*POINTS_PER_DAY);
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_SIMULATE,nStartDate,nEndDate,var,0);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,dFlatRegulationPro);
		app_out_ex_set(4,1);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		memcpy(var,dFlatRegulation,sizeof(double)*POINTS_PER_DAY);
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_SIMULATE,nStartDate,nEndDate,var,0);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,dNegRegulationPro);
		app_out_ex_set(4,2);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		memcpy(var,dNegRegulation,sizeof(double)*POINTS_PER_DAY);
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_SIMULATE,nStartDate,nEndDate,var,0);
		//输出理论典型出力
		app_stat_wind_daily_curve(me, dAverageLoad,dapower,nCount,dPowerMaxall,dPosRegulation,dNegRegulation,dMaxHourlyOutput,dMinHourlyOutput,dAverageHourlyOutput);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,dPosRegulationPro);
		app_out_ex_set(4,0);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		memcpy(var,dPosRegulation,sizeof(double)*POINTS_PER_DAY);
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_THEORETICAL,nStartDate,nEndDate,var,0);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,dFlatRegulationPro);
		app_out_ex_set(4,1);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		memcpy(var,dAverageHourlyOutput,sizeof(double)*POINTS_PER_DAY);
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_THEORETICAL,nStartDate,nEndDate,var,0);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,dNegRegulationPro);
		app_out_ex_set(4,2);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		memcpy(var,dNegRegulation,sizeof(double)*POINTS_PER_DAY);
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_TYPICAL_THEORETICAL,nStartDate,nEndDate,var,0);
	}
	
	free(daLoad);
	free(dapower);
	free(dapowerall);
	free(napowerall);
}

//统计风电预测误差。本软件中统计误差均是预测-实际，正数表示预测高了，需要正备用，否则需要负备用。
void app_out_wind_fpower_stat( app* me)
{
	int i,j,nCount,l,nCount2;
	int nDate,nYear,ntemp;
	unit *pUnit;
	unit_para *pUnitPara,*pUnitParatemp;
	wind_para* pWindPara,*pWindParatemp;
	double *dapower;
	double *dapowerall;
	double *davariation;//里面存各风电场总体variation
	double *davariationall;//里面存所有风电场总体variation
	int* napowerall;
	double var[POINTS_PER_DAY]={0};
	int nBeginYear=0,nEndYear=0;
	double daPos[WPS_DISTRIBUTION_INTERVAL]={0};
	double daNeg[WPS_DISTRIBUTION_INTERVAL]={0};
	double daParameter[6]={0};//里面是6个备用关键参数，分别是正调峰、反调峰以及平均情形下峰荷正备用以及谷荷负备用的值。
	double dPowerMax,dPowerMaxall,dmaxoutputSigma;
	double *daLoad;
	load* pLoad;
	unit_power* pUnitPower,* pUnitFPower;
	double dLoadMax;

	dapower=(double*)malloc(sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
	dapowerall=(double*)malloc(sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
	davariation=(double*)malloc(sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
	davariationall=(double*)malloc(sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
	daLoad=(double*)malloc(sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
	napowerall=(int*)malloc(sizeof(int)*POINTS_PER_DAY*DAYS_PER_YEAR);
	
	nBeginYear=(me->nStartDate/10000);
	if (me->nEndDate==idate_first_day_of_year(me->nEndDate))
	{
		nEndYear=(me->nEndDate/10000)-1;
	}
	else
	{
		nEndYear=(me->nEndDate/10000);
	}
	for (nYear=nBeginYear;nYear<=nEndYear;nYear++)
	{
		dPowerMaxall=0;dmaxoutputSigma=0;
		memset(dapower,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
		memset(dapowerall,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
		memset(davariation,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
		memset(davariationall,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
		memset(daLoad,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
		memset(napowerall,0,sizeof(int)*POINTS_PER_DAY*DAYS_PER_YEAR);

		for(i=0;i<me->aUnits.n;i++)
		{
			memset(dapower,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
			pUnit=(unit*)(me->aUnits.buf[i]);
			nCount=0;dPowerMax=0;
			for (nDate=nYear*10000+101;nDate<nYear*10000+10101;nDate=idate_next_day(nDate))
			{
				if (nDate<me->nStartDate||nDate>=me->nEndDate) continue;
				pUnitPara=app_cur_unit_para(pUnit,nDate);
				if(pUnitPara==NULL) continue;
				if(pUnitPara->nUnitTypeIn!=11) continue;		
				pWindPara=app_cur_wind_para(pUnit,nDate);
				if (pWindPara==NULL) continue;	
				if (dPowerMax<pUnitPara->dPowerMax)
				{
					dPowerMax=pUnitPara->dPowerMax;
				}
				pUnitParatemp=pUnitPara;pWindParatemp=pWindPara;
				pUnitPower=app_cur_unit_power(pUnit,nDate);
				if (pUnitPower==NULL) continue;	
				pUnitFPower=app_cur_unit_fpower(pUnit,nDate);
				if (pUnitFPower==NULL) continue;	
				ntemp=idate_day_idx_of_year(nDate);
				for(j=0;j<POINTS_PER_DAY;j++)	
				{
					dapower[nCount]=pUnitFPower->dPower[j];
					davariation[nCount]=pUnitFPower->dPower[j]-pUnitPower->dPower[j];
					davariationall[(ntemp-1)*POINTS_PER_DAY+j]+=pUnitFPower->dPower[j]-pUnitPower->dPower[j];
					dapowerall[(ntemp-1)*POINTS_PER_DAY+j]+=pUnitFPower->dPower[j];
					napowerall[(ntemp-1)*POINTS_PER_DAY+j]=1;
					nCount++;
				}
			}
			//开始统计风电特性
			if (nCount>0)
			{
				dPowerMaxall+=dPowerMax;
			}
		}
		nCount=0;
		memset(dapower,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
		for(j=0;j<POINTS_PER_DAY*DAYS_PER_YEAR;j++)	
		{
			if (napowerall[j]==1)
			{
				dapower[nCount]=dapowerall[j];
				davariation[nCount]=davariationall[j];
				nCount++;
			}
		}

		dLoadMax=0;
		if (nCount>0)
		{
			memset(daPos,0,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
			memset(daNeg,0,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
			memset(daParameter,0,sizeof(double)*6);
			memset(daLoad,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
			nCount2=0;
			for (nDate=nYear*10000+101;nDate<nYear*10000+10101;nDate=idate_next_day(nDate))
			{
				pLoad=app_cur_load(me,nDate);
				for(j=0;j<POINTS_PER_DAY;j++)	
				{
					daLoad[nCount2+j]=pLoad->dLoad[j];
				}
				for(i=0;i<me->aUnits.n;i++)
				{
					pUnit=(unit*)(me->aUnits.buf[i]);
					pUnitPara=app_cur_unit_para(pUnit,nDate);
					if (pUnitPara!=NULL&&pUnitPara->nUnitTypeIn==11) continue;	
					pUnitPower=app_cur_unit_power(pUnit,nDate);
					if (pUnitPower==NULL) continue;	
					for(j=0;j<POINTS_PER_DAY;j++)	
					{
						if (pUnitPower->dPower[j]<0)
						{
							daLoad[nCount2+j]-=pUnitPower->dPower[j];
						}
					}
				}
				nCount2+=POINTS_PER_DAY;
			}
			for (l=0;l<nCount;l++)
			{
				if (dLoadMax<daLoad[l])
					dLoadMax=daLoad[l];
			}
			//输出风电场总体统计参数
			app_stat_wind_variation(dapower,davariation,daLoad,nCount,dPowerMaxall,daPos,daNeg,daParameter);
			//输出计算结果
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,0);
			app_out_ex_set(4,0);
			app_out_ex_set(5,0);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			var[0]=daParameter[0];var[2]=daParameter[1];
			var[1]=daParameter[0]/dLoadMax;var[3]=daParameter[1]/dLoadMax;
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE_TYPICAL,nYear*10000+101,nYear*10000+10101,var,0);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,0);
			app_out_ex_set(4,3);
			app_out_ex_set(5,0);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			var[0]=daParameter[2];var[2]=daParameter[3];
			var[1]=daParameter[2]/dLoadMax;var[3]=daParameter[3]/dLoadMax;
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE_TYPICAL,nYear*10000+101,nYear*10000+10101,var,0);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,0);
			app_out_ex_set(4,2);
			app_out_ex_set(5,0);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			var[0]=daParameter[4];var[2]=daParameter[5];
			var[1]=daParameter[4]/dLoadMax;var[3]=daParameter[5]/dLoadMax;
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE_TYPICAL,nYear*10000+101,nYear*10000+10101,var,0);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,0);
			app_out_ex_set(4,10);
			app_out_ex_set(5,0);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			memcpy(var,daPos,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE,nYear*10000+101,nYear*10000+10101,var,0);
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,dPowerMaxall);
			app_out_ex_set(3,0);
			app_out_ex_set(4,11);
			app_out_ex_set(5,0);
			memset(var,0,sizeof(double)*POINTS_PER_DAY);
			memcpy(var,daNeg,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE,nYear*10000+101,nYear*10000+10101,var,0);
		}
	}
	free(dapower);
	free(dapowerall);
	free(davariation);
	free(davariationall);
	free(daLoad);
	free(napowerall);
}

//统计风电预测误差。本软件中统计误差均是预测-实际，正数表示预测高了，需要正备用，否则需要负备用。
void app_out_wind_fpower_stat_month( app* me, mtx* pmWindPower[WPS_SIMULATION_TIMES],\
																mtx* pmWindfPower[WPS_SIMULATION_TIMES],int nStartDate, int nEndDate)
{
	int i,j,nCount,k,l,nWindFarm,nCount2;
	int nDate,nmonth;
	unit *pUnit;
	unit_para *pUnitPara;
	wind_para* pWindPara;
	double *dapower;
	double *dapowerall;
	double *davariation;//里面存各风电场总体variation
	double *davariationall;//里面存所有风电场总体variation
	int* napowerall;
	double var[POINTS_PER_DAY]={0};
	int nBeginYear=0,nEndYear=0;
	double daPos[WPS_DISTRIBUTION_INTERVAL]={0};
	double daNeg[WPS_DISTRIBUTION_INTERVAL]={0};
	double daParameter[6]={0};//里面是6个备用关键参数，分别是正调峰、反调峰以及平均情形下峰荷正备用以及谷荷负备用的值。
	double dPowerMax,dPowerMaxall,dmaxoutputSigma;
	double *daLoad;
	load* pLoad;
	unit_power* pUnitPower;
	double dLoadMax;

	dapower=(double*)malloc(sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	dapowerall=(double*)malloc(sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	davariation=(double*)malloc(sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	davariationall=(double*)malloc(sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	daLoad=(double*)malloc(sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	napowerall=(int*)malloc(sizeof(int)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);

	//由于CDF和PDF都需要通过var输出，所以这里规定其长度不能大于var
	if (WPS_DISTRIBUTION_INTERVAL>POINTS_PER_DAY-4)
	{
		app_err(WPS_EXIT_PARAERROR,"输入参数有误");
	}

	nmonth= idate_month_idx_of_year(nStartDate);
	nDate=nStartDate;

	dPowerMaxall=0;dmaxoutputSigma=0;
	memset(dapower,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	memset(dapowerall,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	memset(davariation,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	memset(davariationall,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	memset(daLoad,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	memset(napowerall,0,sizeof(int)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
	nWindFarm=0;
	for(i=0;i<me->aUnits.n;i++)
	{
		memset(dapower,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
		pUnit=(unit*)(me->aUnits.buf[i]);
		pUnitPara=app_cur_unit_para(pUnit,nDate);
		if (pUnitPara==NULL) continue;
		if(pUnitPara->nUnitTypeIn!=11)continue;		
		pWindPara=app_cur_wind_para(pUnit,nDate);
		if (pWindPara==NULL) continue;	
		nCount=0;
		dPowerMax=pUnitPara->dPowerMax;
		nWindFarm++;
		for (l=0;l<WPS_SIMULATION_TIMES;l++)
		{
			for (k=1;k<= pmWindPower[l]->n;k++)
			{
				dapower[nCount]=mtx_get_element(pmWindfPower[l],nWindFarm,k);
				dapowerall[nCount]+=mtx_get_element(pmWindfPower[l],nWindFarm,k);
				davariation[nCount]=mtx_get_element(pmWindfPower[l],nWindFarm,k)-mtx_get_element(pmWindPower[l],nWindFarm,k);
				davariationall[nCount]+=mtx_get_element(pmWindfPower[l],nWindFarm,k)-mtx_get_element(pmWindPower[l],nWindFarm,k);
				napowerall[nCount]=1;
				nCount++;
			}
		}
		if (nCount>0)
		{
			dPowerMaxall+=dPowerMax;
		}
	}
	memset(dapower,0,sizeof(double)*POINTS_PER_DAY*DAYS_PER_YEAR);
	for(j=0;j<nCount;j++)	
	{
		if (napowerall[j]==1)
		{
			dapower[j]=dapowerall[j];
			davariation[j]=davariationall[j];
		}
	}
	dLoadMax=0;
	if (nCount>0)
	{
		//准备各种参数
		memset(daPos,0,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
		memset(daNeg,0,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
		memset(daParameter,0,sizeof(double)*6);
		memset(daLoad,0,sizeof(double)*POINTS_PER_DAY*31*WPS_SIMULATION_TIMES);
		nCount2=0;
		for (l=0;l<WPS_SIMULATION_TIMES;l++)
		{
			for (nDate=nStartDate;nDate<nEndDate;nDate=idate_next_day(nDate))
			{
				pLoad=app_cur_load(me,nDate);
				for(j=0;j<POINTS_PER_DAY;j++)	
				{
					daLoad[nCount2+j]=pLoad->dLoad[j];
				}
				for(i=0;i<me->aUnits.n;i++)
				{
					pUnit=(unit*)(me->aUnits.buf[i]);
					pUnitPara=app_cur_unit_para(pUnit,nDate);
					if (pUnitPara!=NULL&&pUnitPara->nUnitTypeIn==11) continue;	
					pUnitPower=app_cur_unit_power(pUnit,nDate);
					if (pUnitPower==NULL) continue;	
					for(j=0;j<POINTS_PER_DAY;j++)	
					{
						if (pUnitPower->dPower[j]<0)
						{
							daLoad[nCount2+j]-=pUnitPower->dPower[j];
						}
					}
				}
				nCount2+=POINTS_PER_DAY;
			}
		}
		for (l=0;l<nCount;l++)
		{
			if (dLoadMax<daLoad[l])
				dLoadMax=daLoad[l];
		}
		//输出风电场总体统计参数
		app_stat_wind_variation(dapower,davariation,daLoad,nCount,dPowerMaxall,daPos,daNeg,daParameter);
		//输出计算结果
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,0);
		app_out_ex_set(4,0);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		var[0]=daParameter[0];var[2]=daParameter[1];
		var[1]=daParameter[0]/dLoadMax;var[3]=daParameter[1]/dLoadMax;
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE_TYPICAL,nStartDate,nEndDate,var,0);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,0);
		app_out_ex_set(4,3);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		var[0]=daParameter[2];var[2]=daParameter[3];
		var[1]=daParameter[2]/dLoadMax;var[3]=daParameter[3]/dLoadMax;
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE_TYPICAL,nStartDate,nEndDate,var,0);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,0);
		app_out_ex_set(4,2);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		var[0]=daParameter[4];var[2]=daParameter[5];
		var[1]=daParameter[4]/dLoadMax;var[3]=daParameter[5]/dLoadMax;
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE_TYPICAL,nStartDate,nEndDate,var,0);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,0);
		app_out_ex_set(4,10);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		memcpy(var,daPos,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE,nStartDate,nEndDate,var,0);
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,dPowerMaxall);
		app_out_ex_set(3,0);
		app_out_ex_set(4,11);
		app_out_ex_set(5,nmonth);
		memset(var,0,sizeof(double)*POINTS_PER_DAY);
		memcpy(var,daNeg,sizeof(double)*WPS_DISTRIBUTION_INTERVAL);
		app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_ALL_RESERVE,nStartDate,nEndDate,var,0);
	}

	free(dapower);
	free(dapowerall);
	free(davariation);
	free(davariationall);
	free(daLoad);
	free(napowerall);
}
//统计风电预测误差对应正负备用，daParameter里面分别是正调峰峰荷正备用，谷荷负备用，平均出力下峰荷正备用，谷荷负备用，反调峰峰荷正备用，谷荷负备用。
void	app_stat_wind_variation(double *dapower,double *davariation,double *daLoad,int nCount, double dPowerMaxall,\
								double *daPos,double *daNeg,double *daParameter)
{
	int i,j,l,nCount2,nDay;
	double * daVariationtemp,*dapowertemp;
	double dMaxPos,dMaxNeg,dmeanPos,dmeanNeg;
	double dminoutput,dmaxoutput,dmeanoutput;
	double dDaily[POINTS_PER_DAY]={0};
	double daPDF[WPS_DISTRIBUTION_INTERVAL]={0};
	double daCDF[WPS_DISTRIBUTION_INTERVAL]={0};
	int nPeakPoint,nValleyPoint;
	//首先统计各预测水平下风电需求的正负备用值
	daVariationtemp=(double*)malloc(sizeof(double)*nCount);
	dapowertemp=(double*)malloc(sizeof(double)*nCount);
	memset(daVariationtemp,0,sizeof(double)*nCount);
	for (i=0;i<WPS_DISTRIBUTION_INTERVAL;i++)
	{
		nCount2=0;
		for (l=0;l<nCount;l++)
		{
			if ((dapower[l]/(dPowerMaxall+TINY)<((double)(i+1))/WPS_DISTRIBUTION_INTERVAL)&&(dapower[l]/(dPowerMaxall+TINY)>=((double)(i))/WPS_DISTRIBUTION_INTERVAL))
			{
				daVariationtemp[nCount2]=davariation[l];
				nCount2++;
				//printf("预测值%f,预测误差%f\n",dapower[l],davariation[l]);
			}
		}
		app_stat_wind_reserve(daVariationtemp,nCount2,&dMaxPos,&dMaxNeg,&dmeanPos,&dmeanNeg);
		daPos[i]=dmeanPos;
		daNeg[i]=dmeanNeg;
		//printf("%f\n",dmeanPos);
	}

	//统计峰荷时风电出力，确定峰荷时的最大最小
	nDay=(int)nCount/POINTS_PER_DAY;
	nCount2=0;
	memset(daVariationtemp,0,sizeof(double)*nCount);
	memset(dapowertemp,0,sizeof(double)*nCount);
	for (i=0;i<nDay;i++)
	{
		for(j=0;j<POINTS_PER_DAY;j++)
		{
			dDaily[j]=daLoad[nCount2+j];
		}
		nPeakPoint=(int)app_stat_load_peak_valley(dDaily,POINTS_PER_DAY,11);
		dapowertemp[i]=dapower[nPeakPoint+nCount2];
		daVariationtemp[i]=davariation[nPeakPoint+nCount2];
		nCount2+=POINTS_PER_DAY;
	}
	app_stat_wind_power(dapowertemp,nDay,dPowerMaxall,daPDF,daCDF,&dminoutput,&dmaxoutput,&dmeanoutput,WPS_MIN_REGULATION_CONFIDENCE,WPS_MAX_REGULATION_CONFIDENCE);
	app_stat_wind_reserve(daVariationtemp,nDay,&dMaxPos,&dMaxNeg,&dmeanPos,&dmeanNeg);
	daParameter[0]=daPos[(int)(dmaxoutput/dPowerMaxall*WPS_DISTRIBUTION_INTERVAL)];//统计正调峰峰荷正备用容量
	daParameter[4]=daPos[(int)(dminoutput/dPowerMaxall*WPS_DISTRIBUTION_INTERVAL)];//统计反调峰峰荷正备用容量
	daParameter[2]=dmeanPos;//统计峰荷平均正备用容量
	//统计谷荷是的风电出力，确定谷荷时最大最小
	nCount2=0;
	memset(daVariationtemp,0,sizeof(double)*nCount);
	memset(dapowertemp,0,sizeof(double)*nCount);
	for (i=0;i<nDay;i++)
	{
		for(j=0;j<POINTS_PER_DAY;j++)
		{
			dDaily[j]=daLoad[nCount2+j];
		}
		nValleyPoint=(int)app_stat_load_peak_valley(dDaily,POINTS_PER_DAY,12);
		dapowertemp[i]=dapower[nValleyPoint+nCount2];
		daVariationtemp[i]=davariation[nValleyPoint+nCount2];
		nCount2+=POINTS_PER_DAY;
	}
	app_stat_wind_power(dapowertemp,nDay,dPowerMaxall,daPDF,daCDF,&dminoutput,&dmaxoutput,&dmeanoutput,WPS_MIN_REGULATION_CONFIDENCE,WPS_MAX_REGULATION_CONFIDENCE);
	app_stat_wind_reserve(daVariationtemp,nDay,&dMaxPos,&dMaxNeg,&dmeanPos,&dmeanNeg);
	daParameter[1]=daNeg[(int)(dminoutput/dPowerMaxall*WPS_DISTRIBUTION_INTERVAL)];//统计峰荷正调峰谷荷负备用容量
	daParameter[5]=daNeg[(int)(dmaxoutput/dPowerMaxall*WPS_DISTRIBUTION_INTERVAL)];//统计峰荷反调峰谷荷负备用容量
	daParameter[3]=dmeanNeg;//统计谷荷平均负备用容量
	for (i=0;i<WPS_DISTRIBUTION_INTERVAL;i++)
	{
		daPos[i]=daPos[i]/dPowerMaxall;
		daNeg[i]=daNeg[i]/dPowerMaxall;
	}

	free(daVariationtemp);
	free(dapowertemp);
}

void	app_stat_wind_reserve(double *davariation,int nCount, double* dMaxPos,double *dMaxNeg,double *dmeanPos, double *dmeanNeg)
{
	int i,j;
	double *davariationtemp;
	double dtemp;

	if (nCount==0)
	{
		*dMaxPos=0;*dMaxNeg=0;*dmeanPos=0; *dmeanNeg=0;
	}
	else
	{
		davariationtemp=(double*)malloc(sizeof(double)*nCount);
		memset(davariationtemp,0,sizeof(double)*nCount);
		//对出力进行冒泡排序，从小到大
		memcpy(davariationtemp,davariation,sizeof(double)*nCount);
		for(i=0;i<nCount-1;i++)	
		{
			for(j=0;j<nCount-1-i;j++)	
			{
				if (davariationtemp[j]>davariationtemp[j+1])
				{
					dtemp=davariationtemp[j];
					davariationtemp[j]=davariationtemp[j+1];
					davariationtemp[j+1]=dtemp;
				}
			}
		}
		*dMaxNeg=-davariationtemp[util_double2int((1-WPS_MIN_RESERVE_CONFIDENCE)*nCount)];
		j=util_double2int(WPS_MAX_RESERVE_CONFIDENCE*nCount);	
		if (j==nCount) j--;
		*dMaxPos=davariationtemp[j];
		*dmeanPos=0; *dmeanNeg=0;
		for(i=0;i<nCount;i++)	
		{
			if (davariationtemp[i]>0)
				*dmeanPos+=davariationtemp[i]/nCount;
			else 
				*dmeanNeg-=davariationtemp[i]/nCount;
		}
		free(davariationtemp);
	}
}

double app_stat_load_peak_valley(double *daDailyLoad,int nCount,int nType)
{
	int j;
	int nPeakPoint=0,nValleyPoint=0;
	double dPeakLoad=-INFINITY,dValleyLoad=INFINITY;
	for(j=0;j<nCount;j++)	
	{
		if (dPeakLoad<daDailyLoad[j])
		{
			dPeakLoad=daDailyLoad[j];
			nPeakPoint=j;
		}
		else if (dValleyLoad>daDailyLoad[j])
		{
			dValleyLoad=daDailyLoad[j];
			nValleyPoint=j;
		}
	}
	if (nType==0)
		return daDailyLoad[nPeakPoint]-daDailyLoad[nValleyPoint];
	else if (nType==1)
		return daDailyLoad[nPeakPoint];
	else if (nType==2)
		return daDailyLoad[nValleyPoint];
	else if (nType==11)
		return (double)nPeakPoint;
	else if (nType==12)
		return (double)nValleyPoint;
	else
		return 0;
}

double app_wind_variation(double *daWind,int nCount, double dAverage)
{
	int j;
	double dVariation=0;
	for(j=0;j<nCount;j++)	
	{
		dVariation+=fabs(daWind[j]-dAverage);
	}
	return dVariation;
}

void app_stat_wind_typical_curve(app* me, double *daLoad,double *dapower,int nCount,\
	double *dPosRegulation,double *dNegRegulation,double * dFlatRegulation,\
	double *dPosRegulationPro,double *dNegRegulationPro,double *dFlatRegulationPro)
{
	int i,j,nCountNow;
	double daDailyLoad[POINTS_PER_DAY];
	double daDailyLoadWithWind[POINTS_PER_DAY];
	double daWindPower[POINTS_PER_DAY];
	int nDay;
	double dPeakValley,dPeakValleyWithWind,dRatioMax,dRatioMin,dRatiotemp,dWindMean=0,dVariation,dVariationMin=INFINITY;
	memset(dPosRegulation,0,sizeof(double)*POINTS_PER_DAY);
	memset(dNegRegulation,0,sizeof(double)*POINTS_PER_DAY);
	memset(dFlatRegulation,0,sizeof(double)*POINTS_PER_DAY);
	(*dPosRegulationPro)=0;(*dNegRegulationPro)=0;(*dFlatRegulationPro)=0;

	nDay=(int)nCount/POINTS_PER_DAY;
	for(j=0;j<nCount;j++)	
	{
		dWindMean+=dapower[j]/nCount;
	}
	//对日进行循环，分别判断正调峰、逆调峰典型曲线
	nCountNow=0;
	dRatioMax=1;
	dRatioMin=1;
	for (i=0;i<nDay;i++)
	{
		for(j=0;j<POINTS_PER_DAY;j++)	
		{
			daDailyLoad[j]=daLoad[nCountNow+j];
			daDailyLoadWithWind[j]=daLoad[nCountNow+j]-dapower[nCountNow+j];
			daWindPower[j]=dapower[nCountNow+j];
		}
		dPeakValley=app_stat_load_peak_valley(daDailyLoad,POINTS_PER_DAY,0);
		dPeakValleyWithWind=app_stat_load_peak_valley(daDailyLoadWithWind,POINTS_PER_DAY,0);
		dVariation=app_wind_variation(daWindPower,POINTS_PER_DAY, dWindMean);
		dRatiotemp=dPeakValleyWithWind/(dPeakValley+TINY);
		if (dRatiotemp>1+me->dPeakValleyRegulationLimit)
			(*dNegRegulationPro)+=1.0/nDay;
		else if (dRatiotemp<1-me->dPeakValleyRegulationLimit)
			(*dPosRegulationPro)+=1.0/nDay;
		else
			(*dFlatRegulationPro)+=1.0/nDay;

		if (dRatiotemp>dRatioMax)
		{
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				dNegRegulation[j]=dapower[nCountNow+j];
			}
			dRatioMax=dRatiotemp;
		}
		else if (dRatiotemp<dRatioMin)
		{
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				dPosRegulation[j]=dapower[nCountNow+j];
			}
			dRatioMin=dRatiotemp;
		}

		if (dVariation<dVariationMin)
		{
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				dFlatRegulation[j]=dapower[nCountNow+j];
			}
			dVariationMin=dVariation;
		}
		nCountNow+=POINTS_PER_DAY;
	}
}

//设定生成日典型曲线时当负荷接近最高负荷或最低负荷一定裕度下，风电直接置最大或最小出力。
#define WPS_WIND_ALL_TYPICAL_SIMULATE_LIMIT 0.1
//设置理论正调峰反调峰出力曲线
void app_stat_wind_daily_curve(app* me, double *dAverageLoad,double *dapower,int nCount,double dPowerMaxall,\
	double *dPosRegulation,double *dNegRegulation,double *dMaxOutput,\
	double *dMinOutput,double *dAverageOutput)
{
	int j,i,nDays;
	double dMax,dMin,dPeakValley;
	double daPDF[WPS_DISTRIBUTION_INTERVAL]={0};
	double daCDF[WPS_DISTRIBUTION_INTERVAL]={0};
	double dmaxoutput,dminoutput,dmeanoutput;
	double *daHourlyPower;

	nDays=(int)nCount/POINTS_PER_DAY;
	
	app_stat_wind_power(dapower,nCount, dPowerMaxall,daPDF,daCDF,&dminoutput,&dmaxoutput,&dmeanoutput,WPS_MIN_REGULATION_CONFIDENCE,WPS_MAX_REGULATION_CONFIDENCE);
	dMax=app_stat_load_peak_valley(dAverageLoad,POINTS_PER_DAY,1);
	dMin=app_stat_load_peak_valley(dAverageLoad,POINTS_PER_DAY,2);
	dPeakValley=app_stat_load_peak_valley(dAverageLoad,POINTS_PER_DAY,0);
	for(j=0;j<POINTS_PER_DAY;j++)	
	{
		if (dAverageLoad[j]-dMin<WPS_WIND_ALL_TYPICAL_SIMULATE_LIMIT*dPeakValley)
		{
			dNegRegulation[j]=dmaxoutput;dPosRegulation[j]=dminoutput;
		}
		else if (dMax-dAverageLoad[j]<WPS_WIND_ALL_TYPICAL_SIMULATE_LIMIT*dPeakValley)
		{
			dNegRegulation[j]=dminoutput;dPosRegulation[j]=dmaxoutput;
		}
		else
		{
			dNegRegulation[j]=dmaxoutput-(dmaxoutput-dminoutput)*(dAverageLoad[j]-dMin)/(dPeakValley+TINY);
			dPosRegulation[j]=dminoutput+(dmaxoutput-dminoutput)*(dAverageLoad[j]-dMin)/(dPeakValley+TINY);
		}
	}
	daHourlyPower=(double*)malloc(sizeof(double)*nDays);
	memset(daHourlyPower,0,sizeof(double)*nDays);
	for(j=0;j<POINTS_PER_DAY;j++)	
	{
		for(i=0;i<nDays;i++)	
		{
			daHourlyPower[i]=dapower[i*POINTS_PER_DAY+j];
		}
		app_stat_wind_power(daHourlyPower,nDays, dPowerMaxall,daPDF,daCDF,&dminoutput,&dmaxoutput,&dmeanoutput,WPS_MIN_REGULATION_CONFIDENCE,WPS_MAX_REGULATION_CONFIDENCE);
		dMaxOutput[j]=dmaxoutput;
		dMinOutput[j]=dminoutput;
		dAverageOutput[j]=dmeanoutput;
		memset(daHourlyPower,0,sizeof(double)*nDays);
	}
	free(daHourlyPower);
}

void app_out_wind_power( app* me ) 
{
	int i,j;
	int nDate;
	unit *pUnit;
	unit_para *pUnitPara;
	wind_para* pWindPara;
	unit_power* pUnitPower;
	double var[POINTS_PER_DAY]={0};
	double dWindPowerOutput=0;
	memset(var,0,sizeof(double)*POINTS_PER_DAY);
	for (nDate=me->nStartDate;nDate<me->nEndDate;nDate=idate_next_day(nDate))
	{
		for(i=0;i<me->aUnits.n;i++)
		{
			pUnit=(unit*)(me->aUnits.buf[i]);
			pUnitPara=app_cur_unit_para(pUnit,nDate);
			if(pUnitPara==NULL) continue;
			if(pUnitPara->nUnitTypeIn!=11) continue;		
			pWindPara=app_cur_wind_para(pUnit,nDate);
			if (pWindPara==NULL) continue;	
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,pUnitPara->dPowerMax);
			app_out_ex_set(3,pUnitPara->nNodeId);
			app_out_ex_set(4,pWindPara->nWindZone);

			pUnitPower=app_cur_unit_power(pUnit,nDate);
			if (pUnitPower==NULL) continue;	
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				var[j]=pUnitPower->dPower[j];
			}
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_SMLT,nDate,idate_next_day(nDate),var,0);
		}
	}
}

void app_out_wind_fpower( app* me ) 
{
	int i,j;
	int nDate;
	unit *pUnit;
	unit_para *pUnitPara;
	wind_para* pWindPara;
	unit_power* pUnitPower;
	double var[POINTS_PER_DAY]={0};
	double dWindPowerOutput=0;
	memset(var,0,sizeof(double)*POINTS_PER_DAY);
	for (nDate=me->nStartDate;nDate<me->nEndDate;nDate=idate_next_day(nDate))
	{
		for(i=0;i<me->aUnits.n;i++)
		{
			pUnit=(unit*)(me->aUnits.buf[i]);
			pUnitPara=app_cur_unit_para(pUnit,nDate);
			if(pUnitPara==NULL) continue;
			if(pUnitPara->nUnitTypeIn!=11) continue;		
			pWindPara=app_cur_wind_para(pUnit,nDate);
			if (pWindPara==NULL) continue;	
			app_out_ex_clear();
			app_out_ex_set(1,me->nTimes);
			app_out_ex_set(2,pUnitPara->dPowerMax);
			app_out_ex_set(3,pUnitPara->nNodeId);
			app_out_ex_set(4,pWindPara->nWindZone);

			pUnitPower=app_cur_unit_fpower(pUnit,nDate);
			if (pUnitPower==NULL) continue;	
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				var[j]=pUnitPower->dPower[j];
			}
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_FCST,nDate,idate_next_day(nDate),var,0);
		}
	}
}

void app_out_wind_ferror( app* me ) 
{
	int i,j;
	int nDate;
	unit *pUnit;
	unit_para *pUnitPara;
	wind_para* pWindPara;
	unit_power* pUnitPower;
	unit_power* pfUnitPower;
	double var[POINTS_PER_DAY]={0};
	double dt=0;
	double cc=0;
	double dWindPowerOutput=0;
	memset(var,0,sizeof(double)*POINTS_PER_DAY);
	for(i=0;i<me->aUnits.n;i++)
	{
		pUnit=(unit*)(me->aUnits.buf[i]);
		pUnitPara=app_cur_unit_para(pUnit,nDate=me->nStartDate);
		if(pUnitPara==NULL) continue;
		if(pUnitPara->nUnitTypeIn!=11) continue;		
		pWindPara=app_cur_wind_para(pUnit,nDate);
		if (pWindPara==NULL) continue;	
		app_out_ex_clear();
		app_out_ex_set(1,me->nTimes);
		app_out_ex_set(2,pUnitPara->dPowerMax);
		app_out_ex_set(3,pUnitPara->nNodeId);
		app_out_ex_set(4,pWindPara->nWindZone);
		for (nDate;nDate<me->nEndDate;nDate=idate_next_day(nDate))
		{
			pUnitPower=app_cur_unit_power(pUnit,nDate);
			pfUnitPower=app_cur_unit_fpower(pUnit,nDate);
			if (pUnitPower==NULL||pfUnitPower==NULL) continue;
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				var[j]=pfUnitPower->dPower[j]-pUnitPower->dPower[j];
			}
			app_out_array(me,me->nPJID,pUnit->nId,RT_WIND_FERR,nDate,idate_next_day(nDate),var,0);
		}
	}
}