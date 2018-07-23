#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <time.h>
#include <math.h>
#include <stdarg.h>
#include <LIMITS.H>
#include "..\..\BLPK\blp.h"
#include <time.h>
#include "globle.h"
#include "gConstant.h"
#include "func.h"
#pragma warning(disable: 4996)

/*该C文件是进行风电场运行模拟的程序代码，程序中利用风区信息me->aWindZone生成风速序列，考虑风速的Weibull分布参数，
风速季节特性与日特性，风速自相关特性参数以及风电场的相关性。并根据风电场所属风区信息，根据风电场装机容量、风机可靠性参数以及风机尾流效应
生成各风电场出力，风电场出力按月生成，每月生成多组“DESD_WIND_SIMULATION_TIMES个”程序寻找总发电量最接近平均值的那组风电出力
风电预测的误差程序还没有编
*/

/*模拟后一方面将风电出力、预测出力灌入输出csv文件，一方面赋给DSED计算中风电场每日的power、fpower*/

#define DESD_WIND_SIMULATION_TIMES 12/*风电场出力模拟次数（每次模拟均做这么多次，选择其中电量最接近平均电量的那组数据作为最终模拟数据）*/
#define PI 3.1415926

/*以下是一些数学运算函数，取自徐士良《常用算法程序集》*/
//伽马函数
double imaths_gamma(double x)
{ 
	int i;
	double y,t,s,u;
	static double a[11]={ 0.0000677106,-0.0003442342,
		0.0015397681,-0.0024467480,0.0109736958,
		-0.0002109075,0.0742379071,0.0815782188,
		0.4118402518,0.4227843370,1.0};
	if (x<=0.0)
	{ printf("err**x<=0!\n"); return(-1.0);}
	y=x;
	if (y<=1.0)
	{ t=1.0/(y*(y+1.0)); y=y+2.0;}
	else if (y<=2.0)
	{ t=1.0/y; y=y+1.0;}
	else if (y<=3.0) t=1.0;
	else
	{ t=1.0;
	while (y>3.0)
	{ y=y-1.0; t=t*y;}
	}
	s=a[0]; u=y-2.0;
	for (i=1; i<=10; i++)
		s=s*u+a[i];
	s=s*t;
	return(s);
}
//不完全伽马函数
double imaths_gammainc(double a,double x)

{ 
	int n;
	double p,q,d,s,s1,p0,q0,p1,q1,qq;
	if ((a<=0.0)||(x<0.0))
	{ if (a<=0.0) printf("err**a<=0!\n");
	if (x<0.0) printf("err**x<0!\n");
	return(-1.0);
	}
	if (x+1.0==1.0) return(0.0);
	if (x>1.0e+35) return(1.0);
	q=log(x); q=a*q; qq=exp(q);
	if (x<1.0+a)
	{ p=a; d=1.0/a; s=d;
	for (n=1; n<=100; n++)
	{ p=1.0+p; d=d*x/p; s=s+d;
	if (fabs(d)<fabs(s)*1.0e-07)
	{ s=s*exp(-x)*qq/imaths_gamma(a);
	return(s);
	}
	}
	}
	else
	{ s=1.0/x; p0=0.0; p1=1.0; q0=1.0; q1=x;
	for (n=1; n<=100; n++)
	{ p0=p1+(n-a)*p0; q0=q1+(n-a)*q0;
	p=x*p0+n*p1; q=x*q0+n*q1;
	if (fabs(q)+1.0!=1.0)
	{ s1=p/q; p1=p; q1=q;
	if (fabs((s1-s)/s1)<1.0e-07)
	{ s=s1*exp(-x)*qq/imaths_gamma(a);
	return(1.0-s);
	}
	s=s1;
	}
	p1=p; q1=q;
	}
	}
	printf("a too large !\n");
	s=1.0-s*exp(-x)*qq/imaths_gamma(a);
	return(s);
}
//生成均值为u，方差为g^2的正态分布的随机数
double imaths_normal(double u,double g,double* r)
{ 
	int i,m;
	double s,w,v,t;
	s=65536.0; w=2053.0; v=13849.0;
	t=0.0;
	for (i=1; i<=12; i++)
	{ *r=(*r)*w+v; m=(int)(*r/s);
	*r=*r-m*s; t=t+(*r)/s;
	}
	t=u+g*(t-6.0);
	return(t);
}

//生成正态分布的随机数的外层函数，该函数里面定义了静态变量作为随机数的种子
static double imaths_gen_normal_r_wind_speed=3;//专门用于生成风速的随机数种子，为了保持每一年风速的一致性，便于多年之间的比较
static double imaths_gen_normal_r=3;//公共随机数种子
double imaths_gen_normal(double u,double g, int *refresh)
{
	if (refresh!=NULL&&(*refresh)!=1)//如果还原标志refresh不为1，则不重置随机数种子
	{
		return(imaths_normal(u,g,&imaths_gen_normal_r_wind_speed));
	}
	else if (refresh!=NULL&&(*refresh)==1)//如果还原标志refresh置1，则重置随机数种子
	{
		imaths_gen_normal_r_wind_speed=3;
		(*refresh)=0;//重置之后就将其改为0，直至下次被它被外部再改为1
		return(imaths_normal(u,g,&imaths_gen_normal_r_wind_speed));
	}
	else//如果refresh是NULL则，直接用公共随机数种子
	{
		return(imaths_normal(u,g,&imaths_gen_normal_r));
	}
}

//生成参数n，p的贝努力分布的随机数（每次试验成功的概率为p，做n次独立试验成功的总次数服从贝努力分布）
int imaths_bernoulli(int n,double p)
{ 
	double t;
	t=imaths_gen_normal(n*p,n*p*(1-p),NULL);
	return (int)t>n?n:(int)t;
}

//对称正定阵的LL分解
int imaths_chol(double a[],int n,double *det)
{ 
	int i,j,k,u,l;
	double d;
	if ((a[0]+1.0==1.0)||(a[0]<0.0))
	{ printf("fail\n"); return(FALSE);}
	a[0]=sqrt(a[0]);
	d=a[0];
	for (i=1; i<=n-1; i++)
	{ u=i*n; a[u]=a[u]/a[0];}
	for (j=1; j<=n-1; j++)
	{ l=j*n+j;
	for (k=0; k<=j-1; k++)
	{ u=j*n+k; a[l]=a[l]-a[u]*a[u];}
	if ((a[l]+1.0==1.0)||(a[l]<0.0))
	{ printf("fail\n"); return(FALSE);}
	a[l]=sqrt(a[l]);
	d=d*a[l];
	for (i=j+1; i<=n-1; i++)
	{ u=i*n+j;
	for (k=0; k<=j-1; k++)
		a[u]=a[u]-a[i*n+k]*a[j*n+k];
	a[u]=a[u]/a[l];
	}
	}
	*det=d*d;
	for (i=0; i<=n-2; i++)
		for (j=i+1; j<=n-1; j++)
			a[i*n+j]=0.0;
	return(TRUE);
}

//对称正定阵的LL分解,输入mtx数据类型
int imaths_chol_mtx(mtx * pm)
{
	double det;
	if (pm->m!=pm->n)
	{
		printf("pm matrix has to be square matrix  !\n");
		return FALSE;
	}
	else
	{
		return imaths_chol(mtx_get_buf_write(pm),pm->m,&det);
	}
}

//计算x的a次幂，x需大于等于零
double imaths_power(double x, double a)
{
	if (x<0)
	{
		printf("x has to larger than 0  !\n");return 0;
	}
	else if (x==0)
	{
		if (a<0)
		{
			printf("when x==0 a has to larger than 0  !\n");
			return 0;
		}
		else if (a==0)
			return 1.0;
		else 
			return 0.0;
	}
	return exp((log(x)*a));
}

//返回双参数Weibull的x处的密度函数值
double imaths_Weibull_f(double x,double c,double k)
{
	if (x<=0)
		return 0.0;
	else
		return (k/c)*(imaths_power(x/c,k-1))*exp(-(imaths_power(x/c,k)));
}

//返回双参数Weibull的x处的概率分布函数值
double imaths_Weibull_F(double x,double c,double k)
{
	if (x<=0)
		return 0.0;
	else
		return 1-exp(-(imaths_power(x/c,k)));
}


/*以下是风电场运行模拟程序*/

/*主函数控制风电场出力按月生成，每月生成多组
“DESD_WIND_SIMULATION_TIMES个”程序寻找总发电量最接近平均值的那组风电出力*/
void app_generate_wind_power(app *me)
{
	int nDate,nDateStart,nDateEnd,i,nmonth;
	int nWindFarm;
	double dWindFarmCap,dWindFarmRateAve,dWindFarmFRateAve;
	double dWindFarmRate[DESD_WIND_SIMULATION_TIMES];
	double dWindFarmFRate[DESD_WIND_SIMULATION_TIMES];
	int nChoice[12];
	mtx* pmWindPower[DESD_WIND_SIMULATION_TIMES];
	mtx* pmWindFPower[DESD_WIND_SIMULATION_TIMES];
	unit* pUnit;
	unit_para* pUnitPara;
	wind_para* pWindPara;
	int nRefresh;

#define DSED_GENERATE_WIND_POWER_OUT_MATRIX
#ifdef DSED_GENERATE_WIND_POWER_OUT_MATRIX
	char  f_caFullFileName[255];	
	FILE *f1,*f2;
	double var[POINTS_PER_DAY];
	int j;
	unit_power* pUnitPower;
	sprintf(f_caFullFileName,"%swindpower.csv",me->caOutPathName);
	f1= fopen(f_caFullFileName,"w+");
	if(f1==NULL)
	{
		printf("无法写入数据文件：'windpower.csv'");
		return;
	}
	sprintf(f_caFullFileName,"%swindpowerall.csv",me->caOutPathName);
	f2= fopen(f_caFullFileName,"w+");
	if(f2==NULL)
	{
		printf("无法写入数据文件：'windpowerall.csv'");
		return;
	}
#endif

	for (i=0;i<12;i++)
	{
		nChoice[i]=-1;
	}

	printf("【提示信息】风电场出力模拟中……\n");
	for (nDate=me->nStartDate;nDate<me->nEndDate;)
	{
		nWindFarm=0;
		dWindFarmCap=0;
		dWindFarmRateAve=0;
		dWindFarmFRateAve=0;
		nDateStart=nDate;
		nDateEnd=idate_next_month(nDateStart);
		if (nDateEnd>me->nEndDate) nDateEnd=me->nEndDate;

		if (nDateStart==idate_first_day_of_year(nDateStart))
			nRefresh=1;
		else
			nRefresh=0;
		for(i=0;i<me->aUnits.n;i++)
		{
			pUnit=(unit*)(me->aUnits.buf[i]);
			pUnitPara=app_cur_unit_para(pUnit,nDate);
			if (pUnitPara==NULL) continue;
			if(pUnitPara->nUnitTypeIn!=11)continue;		
			pWindPara=app_cur_wind_para(pUnit,nDate);
			if (pWindPara==NULL) continue;	
			pWindPara->nNumofTurbine=util_double2int(pUnitPara->dPowerMax);
			nWindFarm++;
			dWindFarmCap+=pUnitPara->dPowerMax;
		}
		if (nWindFarm==0)
		{
			nDate=nDateEnd;
			continue;
		}
		for (i=0;i<DESD_WIND_SIMULATION_TIMES;i++)
		{
			pmWindPower[i]=mtx_create(nWindFarm,idate_days_between(nDateStart,nDateEnd)*POINTS_PER_DAY,1);
			Generate_wind_power_month(nDateStart,nDateEnd,me,pmWindPower[i],&nRefresh);
			dWindFarmRate[i]=mtx_sum_element(pmWindPower[i])/(dWindFarmCap*idate_days_between(nDateStart,nDateEnd)*POINTS_PER_DAY);
			dWindFarmRateAve+=dWindFarmRate[i]/DESD_WIND_SIMULATION_TIMES;

			pmWindFPower[i]=mtx_create(nWindFarm,idate_days_between(nDateStart,nDateEnd)*POINTS_PER_DAY,1);
			Generate_wind_forecasterror_month(nDateStart,nDateEnd,me,pmWindPower[i],pmWindFPower[i],&nRefresh);
			dWindFarmFRate[i]=mtx_sum_element(pmWindFPower[i])/(dWindFarmCap*idate_days_between(nDateStart,nDateEnd)*POINTS_PER_DAY);
			dWindFarmFRateAve+=dWindFarmFRate[i]/DESD_WIND_SIMULATION_TIMES;
		}
		nmonth=idate_month_idx_of_year(nDate);
		nmonth=nmonth-1;
		if (nChoice[nmonth]==-1)
		{
			nChoice[nmonth]=0;
			for (i=0;i<DESD_WIND_SIMULATION_TIMES;i++)
			{
				if (fabs(dWindFarmRate[i]-dWindFarmRateAve)<fabs(dWindFarmRate[nChoice[nmonth]]-dWindFarmRateAve)-TINY)
					nChoice[nmonth]=i;
			}
		}
		//赋给DSED
		Generate_wind_power_SetWindPower(nDateStart,nDateEnd,me,pmWindPower[nChoice[nmonth]]);
		Generate_wind_power_SetWindFPower(nDateStart,nDateEnd,me,pmWindFPower[nChoice[nmonth]]);

#ifdef DSED_GENERATE_WIND_POWER_OUT_MATRIX
		//按列输出各风电场出力
		mtx_fprint(f1,mtx_transpose_self(pmWindPower[nChoice[nmonth]]));
#endif

		for (i=0;i<DESD_WIND_SIMULATION_TIMES;i++)
		{
			mtx_free(pmWindPower[i]);
			mtx_free(pmWindFPower[i]);
		}
		nDate=nDateEnd;
	}

	Generate_wind_power_StatisticsWindPower(me);

#ifdef DSED_GENERATE_WIND_POWER_OUT_MATRIX
	//按日期（行）输出风电场总出力
	for (nDate=me->nStartDate;nDate<me->nEndDate;nDate=idate_next_day(nDate))
	{
		memset(var,0,sizeof(double)*(POINTS_PER_DAY));
		if (nDate%10000==229) continue;
		for(i=0;i<me->aUnits.n;i++)
		{
			pUnit=(unit*)(me->aUnits.buf[i]);
			pUnitPara=app_cur_unit_para(pUnit,nDate);
			if(pUnitPara==NULL) continue;
			if(pUnitPara->nUnitTypeIn!=11) continue;		
			pUnitPower=app_cur_unit_power(pUnit,nDate);
			if (pUnitPower==NULL) continue;
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				var[j]+=pUnitPower->dPower[j];
			}
		}
		fprintf(f2,"%d,",nDate);
		for(j=0;j<POINTS_PER_DAY;j++)	
		{
			fprintf(f2,"%.10g,",var[j]);
		}
		fprintf(f2,"\n");
	}
	fclose(f1);
	fclose(f2);
#endif

	printf("【提示信息】风电场出力模拟完毕\n");
}

/*按月生成风电场出力的函数*/
void Generate_wind_power_month(int nDateStart,int nDateEnd,app *me,mtx* pmWindPower, int *nRefresh)
{
	int nDate;
	int i,j,l,nbas;
	int nmonth,nWindFarm;
	double c,k,theta;
	double dx,dWSPre,dNu;
	double dtemp;
	unit* pUnit;
	unit_para* pUnitPara;
	wind_para* pWindPara;
	wind_zone* pWindzone;
	mtx* pmDW;
	mtx* pmWindSpeedPre;
	mtx* pmL;
	mtx* pmDWCorr;
	mtx* pmWindSpeed;
	int nFarmAvailable;
	int nWindZone=me->aWindZone.n;
	wind_zone_para* pWindZonePara;

	if (nWindZone==0) return;

	//mtx_print(me->pmWindCorr);
	pmL=mtx_create(0,0,0);
	mtx_copy(pmL,me->pmWindCorr);
	imaths_chol_mtx(pmL);
	//mtx_print(pmL);

	//初始化第一时段前一时段的风速（现在暂且认为是平均风速）
	pmWindSpeedPre=mtx_create(nWindZone,1,1);
	for(i=0;i<nWindZone;i++)
	{
		pWindzone=(wind_zone*)(me->aWindZone.buf[i]);
		pWindZonePara=(wind_zone_para*)pa_get(&pWindzone->para,0);
		c=pWindZonePara->dWeibull_c;
		k=pWindZonePara->dWeibull_k;
		theta=pWindZonePara->dtheta;
		mtx_set_element(pmWindSpeedPre,i+1,1,c*imaths_gamma(1/k+1));
	}
	//mtx_print(pmWindSpeedPre);
	//以下过程按日进行
	pmDW=mtx_create(nWindZone,POINTS_PER_DAY,1);
	pmDWCorr=mtx_create(nWindZone,POINTS_PER_DAY,1);
	pmWindSpeed=mtx_create(nWindZone,POINTS_PER_DAY,1);
	nbas=0;
	for (nDate=nDateStart;nDate<nDateEnd;nDate=idate_next_day(nDate))
	{
		//第一步，生成相关系数矩阵为me->pmWindCorr的正态分布序列
		for(i=0;i<me->aWindZone.n;i++)
		{
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				mtx_set_element(pmDW,i+1,j+1,imaths_gen_normal(0,1,nRefresh));
			}
		}
		//mtx_print(pmDW);
		mtx_mul(pmDWCorr,pmL,pmDW);
		//mtx_print(pmDWCorr);
		//第二步，利用pmWindCorr以及随机差分方程生成标准风速序列
		for(i=0;i<nWindZone;i++)
		{
			pWindzone=(wind_zone*)(me->aWindZone.buf[i]);
			pWindZonePara=(wind_zone_para*)pa_get(&pWindzone->para,0);
			c=pWindZonePara->dWeibull_c;
			k=pWindZonePara->dWeibull_k;
			theta=pWindZonePara->dtheta;
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				if (j==0) dWSPre=mtx_get_element(pmWindSpeedPre,i+1,1);
				else dWSPre=mtx_get_element(pmWindSpeed,i+1,j+1-1);
				dNu=Generate_wind_power_Nu(dWSPre,c,k,theta);
				dx=theta*(c*imaths_gamma(1+1/k)-dWSPre)+imaths_power(dNu,0.5)*mtx_get_element(pmDWCorr,i+1,j+1);
				if (dWSPre+dx<0) mtx_set_element(pmWindSpeed,i+1,j+1,0);
				else mtx_set_element(pmWindSpeed,i+1,j+1,dWSPre+dx);
			}
			mtx_set_element(pmWindSpeedPre,i+1,1,mtx_get_element(pmWindSpeed,i+1,POINTS_PER_DAY));
		}
		//mtx_print(pmWindSpeed);
		//mtx_print(pmWindSpeedPre);
		//第三步，根据风速的季节特性与日特性修正风速
		for(i=0;i<nWindZone;i++)
		{
			pWindzone=(wind_zone*)(me->aWindZone.buf[i]);
			nmonth=idate_month_idx_of_year(nDate);
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				dtemp=mtx_get_element(pmWindSpeed,i+1,j+1);
				pWindZonePara=(wind_zone_para*)pa_get(&pWindzone->para,0);
				dtemp*=pWindZonePara->nDayAve[j]*pWindZonePara->nMonthAve[nmonth-1];
				mtx_set_element(pmWindSpeed,i+1,j+1,dtemp);
			}
		}
		//mtx_print(pmWindSpeed);
		//第四步，根据生成的风速序列计算风电场出力
		nWindFarm=0;
		for(i=0;i<me->aUnits.n;i++)
		{
			pUnit=(unit*)(me->aUnits.buf[i]);
			pUnitPara=app_cur_unit_para(pUnit,nDate);
			if (pUnitPara==NULL) continue;
			if(pUnitPara->nUnitTypeIn!=11)continue;		
			pWindPara=app_cur_wind_para(pUnit,nDate);
			if (pWindPara==NULL) continue;		
			nFarmAvailable=imaths_bernoulli(pWindPara->nNumofTurbine,pWindPara->dAvailableRate);
			//nFarmAvailable=(int)(pWindPara->nNumofTurbine*pWindPara->dAvailableRate);
			nWindFarm++;
			for(l=0;l<nWindZone;l++)
			{
				pWindzone=(wind_zone*)(me->aWindZone.buf[l]);
				if(pWindzone->nId==pWindPara->nWindZone)
				{
					for(j=0;j<POINTS_PER_DAY;j++)	
					{
						dtemp=Generate_wind_power_curve(mtx_get_element(pmWindSpeed,l+1,j+1),
							pWindPara->dCutinSP, pWindPara->dRatedSP, pWindPara->dCutoutSP);
						dtemp*=nFarmAvailable;
						dtemp*=pWindPara->dWakeeffect;
						mtx_set_element(pmWindPower,nWindFarm,nbas+j+1,dtemp);
					}
				}
			}
		}	
		nbas+=POINTS_PER_DAY;
	}
	//mtx_print(pmWindPower);
	mtx_free(pmDW);
	mtx_free(pmWindSpeedPre);
	mtx_free(pmL);	
	mtx_free(pmDWCorr);
	mtx_free(pmWindSpeed);
}

/*按月生成风电场出力预测误差*/
void Generate_wind_forecasterror_month(int nDateStart,int nDateEnd,app *me,mtx* pmWindPower,mtx* pmWindFPower,int *nRefresh)
{
	int nDate;
	int i,j,l,nbas;
	int nWindFarm;
	double dMRE,dWindFarmMax,theta;
	double dx,dWFEPre,dNu;
	double dtemp;
	unit* pUnit;
	unit_para* pUnitPara;
	wind_para* pWindPara;
	wind_zone* pWindzone;
	mtx* pmDW;
	mtx* pmWindZoneErrorPre;
	mtx* pmL;
	mtx* pmDWCorr;
	mtx* pmWindZoneError;
	int nWindZone=me->aWindZone.n;
	wind_zone_para* pWindZonePara;
	if (nWindZone==0) return;


	//mtx_print(me->pmWindCorr);
	pmL=mtx_create(0,0,0);
	mtx_copy(pmL,me->pmWindCorr);
	imaths_chol_mtx(pmL);
	//mtx_print(pmL);

	//初始化前一时段的各风区的预测误差（现在暂且认为是0）
	pmWindZoneErrorPre=mtx_create(nWindZone,1,1);
	for(i=0;i<nWindZone;i++)
	{
		mtx_set_element(pmWindZoneErrorPre,i+1,1,0);
	}
	//mtx_print(pmWindSpeedPre);
	//以下过程按日进行
	pmDW=mtx_create(nWindZone,POINTS_PER_DAY,1);
	pmDWCorr=mtx_create(nWindZone,POINTS_PER_DAY,1);
	pmWindZoneError=mtx_create(nWindZone,POINTS_PER_DAY,1);
	nbas=0;
	for (nDate=nDateStart;nDate<nDateEnd;nDate=idate_next_day(nDate))
	{
		//第一步，生成相关系数矩阵为me->pmWindCorr的正态分布序列
		for(i=0;i<me->aWindZone.n;i++)
		{
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				mtx_set_element(pmDW,i+1,j+1,imaths_gen_normal(0,1,nRefresh));
			}
		}
		//mtx_print(pmDW);
		mtx_mul(pmDWCorr,pmL,pmDW);
		//mtx_print(pmDWCorr);
		//第二步，利用pmWindCorr以及随机差分方程生成预测误差标幺值序列
		for(i=0;i<nWindZone;i++)
		{
			pWindzone=(wind_zone*)(me->aWindZone.buf[i]);
			pWindZonePara=(wind_zone_para*)pa_get(&pWindzone->para,0);
			theta=pWindZonePara->dtheta;
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				if (j==0) dWFEPre=mtx_get_element(pmWindZoneErrorPre,i+1,1);
				else dWFEPre=mtx_get_element(pmWindZoneError,i+1,j+1-1);
				dNu=Generate_wind_forecasterror_Nu(theta);
				dx=theta*(-dWFEPre)+imaths_power(theta*PI,0.5)*mtx_get_element(pmDWCorr,i+1,j+1);
				mtx_set_element(pmWindZoneError,i+1,j+1,dWFEPre+dx);
			}
			mtx_set_element(pmWindZoneErrorPre,i+1,1,mtx_get_element(pmWindZoneError,i+1,POINTS_PER_DAY));
		}
		mtx_sum_element(pmWindZoneError);
		//mtx_print(pmWindZoneError);
		//mtx_print(pmWindFErrorPre);

		//mtx_print(pmWindZoneError);
		//第三步，根据各风区预测标幺误差，确定各风电场预测误差
		nWindFarm=0;
		for(i=0;i<me->aUnits.n;i++)
		{
			pUnit=(unit*)(me->aUnits.buf[i]);
			pUnitPara=app_cur_unit_para(pUnit,nDate);
			if (pUnitPara==NULL) continue;
			if(pUnitPara->nUnitTypeIn!=11)continue;		
			pWindPara=app_cur_wind_para(pUnit,nDate);
			if (pWindPara==NULL) continue;		
			nWindFarm++;
			dMRE=pWindPara->dForecastError;
			dWindFarmMax=pUnitPara->dPowerMax;

			if (pUnit->nId==1121)
			{
				//printf("测试");
			}
			for(l=0;l<nWindZone;l++)
			{
				pWindzone=(wind_zone*)(me->aWindZone.buf[l]);
				if(pWindzone->nId==pWindPara->nWindZone)
				{
					for(j=0;j<POINTS_PER_DAY;j++)	
					{
						dtemp=mtx_get_element(pmWindZoneError,l+1,j+1)*dMRE*dWindFarmMax;
						dtemp=dtemp+mtx_get_element(pmWindPower,nWindFarm,nbas+j+1);
						if (dtemp<TINY) dtemp=0;
						else if (dtemp>dWindFarmMax-TINY) dtemp=dWindFarmMax;
						mtx_set_element(pmWindFPower,nWindFarm,nbas+j+1,dtemp);
					}
				}
			}
		}	
		nbas+=POINTS_PER_DAY;
	}
	//mtx_print(pmWindFPower);
	mtx_free(pmDW);
	mtx_free(pmWindZoneErrorPre);
	mtx_free(pmL);	
	mtx_free(pmDWCorr);
	mtx_free(pmWindZoneError);
}

/*风速模拟中随机差分方程中根号下那个函数，应该读作Ν ν：拗 Nu */
double Generate_wind_power_Nu(double x, double c, double k, double theta)
{
	double F,a,b;

	F=imaths_Weibull_F(x,c,k);
	a=2*theta*c*imaths_gamma(1+1/k)*F-2*theta*c*imaths_gammainc(1+1/k,imaths_power(x/c,k))*imaths_gamma(1+1/k);
	b=imaths_Weibull_f(x,c,k);

	return a/(b+TINY);
}

double Generate_wind_forecasterror_Nu(double theta)
{
	return imaths_power(theta*PI,0.5);
}

//风机功率特性曲线（标幺化）
double Generate_wind_power_curve(double dWindSP, double dCutinSP, double dRatedSP, double dCutoutSP)
{
	if (dWindSP>=dCutinSP&&dWindSP<=dRatedSP)
		return (imaths_power(dWindSP,3)-imaths_power(dCutinSP,3))/(imaths_power(dRatedSP,3)-imaths_power(dCutinSP,3));
	else if (dWindSP>=dRatedSP&&dWindSP<=dCutoutSP)
		return 1.0;
	else
		return 0.0;
}

//将模拟得到的机组出力置每个风电场每日的power（数据类型为unit_power）
void Generate_wind_power_SetWindPower(int nDateStart,int nDateEnd,app *me,mtx* pmWindPower)
{
	int nDate;
	int i,j,nbas;
	int nWindFarm;
	unit* pUnit;
	unit_para* pUnitPara;
	wind_para* pWindPara;
	unit_power* pUnitPower;
	nWindFarm=0;
	for(i=0;i<me->aUnits.n;i++)
	{
		pUnit=(unit*)(me->aUnits.buf[i]);
		if (pUnit->nId==1117)
		{
			//printf("测试");
		}
		pUnitPara=app_cur_unit_para(pUnit,nDateStart);
		if(pUnitPara==NULL) continue;
		if(pUnitPara->nUnitTypeIn!=11) continue;		
		pWindPara=app_cur_wind_para(pUnit,nDateStart);
		if (pWindPara==NULL) continue;	
		nWindFarm++;
		nbas=0;
		for (nDate=nDateStart;nDate<nDateEnd;nDate=idate_next_day(nDate))
		{
			pUnitPower=app_cur_unit_power(pUnit,nDate);
			if (pUnitPower==NULL)
			{
				pUnitPower=malloc(sizeof(unit_power));
				pa_add(&pUnit->power,pUnitPower);
				pUnitPower->nDate=nDate;
				pUnitPower->nDate2=idate_next_day(nDate);
				for(j=0;j<POINTS_PER_DAY;j++)	
				{
					pUnitPower->dPower[j]=mtx_get_element(pmWindPower,nWindFarm,nbas+j+1);
				}
			}
			nbas+=POINTS_PER_DAY;
		}
	}
}

//将模拟得到的机组预测出力置每个风电场每日的fpower（数据类型为unit_power）
void Generate_wind_power_SetWindFPower(int nDateStart,int nDateEnd,app *me,mtx* pmWindFPower)
{
	int nDate;
	int i,j,nbas;
	int nWindFarm;
	unit* pUnit;
	unit_para* pUnitPara;
	wind_para* pWindPara;
	unit_power* pUnitPower;
	nWindFarm=0;
	for(i=0;i<me->aUnits.n;i++)
	{
		pUnit=(unit*)(me->aUnits.buf[i]);
		if (pUnit->nId==1117)
		{
			//printf("测试");
		}
		pUnitPara=app_cur_unit_para(pUnit,nDateStart);
		if(pUnitPara==NULL) continue;
		if(pUnitPara->nUnitTypeIn!=11) continue;		
		pWindPara=app_cur_wind_para(pUnit,nDateStart);
		if (pWindPara==NULL) continue;	
		nWindFarm++;
		nbas=0;
		for (nDate=nDateStart;nDate<nDateEnd;nDate=idate_next_day(nDate))
		{
			pUnitPower=malloc(sizeof(unit_power));
			pa_add(&pUnit->fpower,pUnitPower);
			pUnitPower->nDate=nDate;
			pUnitPower->nDate2=idate_next_day(nDate);
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				pUnitPower->dPower[j]=mtx_get_element(pmWindFPower,nWindFarm,nbas+j+1);
			}
			nbas+=POINTS_PER_DAY;
		}
	}
}

//统计年风电场总出力分布
void Generate_wind_power_StatisticsWindPower(app *me)
{
	int nDate,nDateStart,nDateEnd,i,j,nWindFarm,nPoints;
	double dWindFarmCap;
	double dWindPower[POINTS_PER_DAY];
	unit* pUnit;
	unit_para* pUnitPara;
	unit_power* pUnitPower;
	wind_powerPDF* pWindPowerPDF;

#ifdef DSED_GENERATE_WIND_POWER_OUT_MATRIX
	char  f_caFullFileName[255];	
	FILE *f1;
	sprintf(f_caFullFileName,"%swindpowerPDF.csv",me->caOutPathName);
	f1= fopen(f_caFullFileName,"w+");
	if(f1==NULL)
	{
		printf("无法写入数据文件：'windpower.csv'");
		return;
	}
#endif

	for (nDate=me->nStartDate;nDate<me->nEndDate;)
	{
		nWindFarm=0;
		dWindFarmCap=0;
		nDateStart=nDate;
		nDateEnd=idate_first_day_of_year(nDateStart+10000);
		if (nDateEnd>me->nEndDate) nDateEnd=me->nEndDate;
		for(i=0;i<me->aUnits.n;i++)
		{
			pUnit=(unit*)(me->aUnits.buf[i]);
			pUnitPara=app_cur_unit_para(pUnit,nDate);
			if (pUnitPara==NULL) continue;
			if(pUnitPara->nUnitTypeIn!=11)continue;		
			nWindFarm++;
			dWindFarmCap+=pUnitPara->dPowerMax;
		}
		if (nWindFarm==0)
		{
			nDate=nDateEnd;
			continue;
		}

		pWindPowerPDF=malloc(sizeof(wind_powerPDF));
		memset(pWindPowerPDF,0,sizeof(wind_powerPDF));
		pWindPowerPDF->nDate=nDateStart;
		pWindPowerPDF->nDate2=nDateEnd;
		pWindPowerPDF->PowerPDF=mtx_create((int)(dWindFarmCap/WPS_RELIABILITY_STEP)+1,1,0);
		mtx_ini_zeros(pWindPowerPDF->PowerPDF,(int)(dWindFarmCap/WPS_RELIABILITY_STEP)+1,1);
		pa_add(&me->aWindPowerPDF,pWindPowerPDF);

		nPoints=0;
		for (nDate=nDateStart;nDate<nDateEnd;nDate=idate_next_day(nDate))
		{
			if (nDate==20200128)
			{
				//printf("测试");
			}
			memset(dWindPower,0,sizeof(double)*POINTS_PER_DAY);
			for(i=0;i<me->aUnits.n;i++)
			{
				pUnit=(unit*)(me->aUnits.buf[i]);
				pUnitPara=app_cur_unit_para(pUnit,nDateStart);
				if(pUnitPara==NULL) continue;
				if(pUnitPara->nUnitTypeIn!=11) continue;		
				pUnitPower=app_cur_unit_power(pUnit,nDate);
				if (pUnitPower!=NULL)
				{
					for(j=0;j<POINTS_PER_DAY;j++)	
					{
						dWindPower[j]+=pUnitPower->dPower[j];
					}
				}
			}
			for(j=0;j<POINTS_PER_DAY;j++)	
			{
				if ((int)dWindPower[j]<=(int)dWindFarmCap)
				{
					nPoints++;
					mtx_set_element(pWindPowerPDF->PowerPDF,(int)(dWindPower[j]/WPS_RELIABILITY_STEP+0.4999999999999)+1,1,\
						mtx_get_element(pWindPowerPDF->PowerPDF,(int)(dWindPower[j]/WPS_RELIABILITY_STEP+0.4999999999999)+1,1)+1);
				}
			}
		}
		mtx_sum_element(pWindPowerPDF->PowerPDF);
		mtx_mul_number_self(1.0/nPoints,pWindPowerPDF->PowerPDF);
		mtx_sum_element(pWindPowerPDF->PowerPDF);

#ifdef DSED_GENERATE_WIND_POWER_OUT_MATRIX
		fprintf(f1,"*********%d年**********\n",nDateStart/10000);
		mtx_fprint(f1,pWindPowerPDF->PowerPDF);	
#endif	
		nDate=nDateEnd;
	}
#ifdef DSED_GENERATE_WIND_POWER_OUT_MATRIX
	fclose(f1);
#endif
}

/*计算某一时刻风电容量可信度，其中z是没有包含风电出力情况下机组可调出力概率，
dLOLP为考虑风电之后的LOLP，dFor是等效火电的强迫停运率，nLoad是原始负荷。 dmax是迭代开始可信容量最大值，一般取风电场总容量。
*/
double calc_WindCC_on_LOLPorEENS(mtx* z,double dLOLPorEENS, double dFor,int nLoad, double dmax, int nIsEENS)
{
	double dCapamin=0;
	double dCapamax= dmax;
	double dCapamid;
	double dLE=dLOLPorEENS;
	double dLEmid=0,dLEmax,dLEmin;


	dLEmax=calc_LOLPorEENS(z,dCapamin,dFor,nLoad,nIsEENS);
	if (dLEmax<dLE)
	{
		return dCapamin;
	}
	dLEmin=calc_LOLPorEENS(z,dCapamax,dFor,nLoad,nIsEENS);
	if (dLEmin>dLE)
	{
		return dCapamax;
	}
	dCapamid=(dCapamin+dCapamax)/2;

	dLEmid=calc_LOLPorEENS(z,dCapamid,dFor,nLoad,nIsEENS);
	while (fabs(dLE-dLEmid)/dLE>0.001&&dLEmid!=dLEmax&&dLEmid!=dLEmin)
	{
		if (dLEmid>dLE)
		{
			dCapamin=dCapamid;
			dLEmax=dLEmid;
		}
		else
		{
			dCapamax=dCapamid;
			dLEmin=dLEmid;
		}
		dCapamid=(dCapamin+dCapamax)/2;
		dLEmid=calc_LOLPorEENS(z,dCapamid,dFor,nLoad,nIsEENS);
	}
	return dCapamid;
}

double calc_LOLPorEENS(mtx* z,double dCapa, double dFor,int nLoad,int nIsEENS)
{
	double *data;
	double dt,dt2;
	double dFOR=dFor;
	double d1_FOR=1-dFOR;
	double dLE=0;
	int nSizeOriginal,nSizeplus,m;
	int nUnitCapa=(int)(dCapa/WPS_RELIABILITY_STEP+0.4999999999999);

	nSizeOriginal=z->m;
	nSizeplus=nSizeOriginal+nUnitCapa;
	data=malloc(sizeof(double)*nSizeplus);
	/*计算移位后不和原始序列重叠部分*/
	for(m=nSizeOriginal;m<nSizeplus;m++)
	{
		if (m>=nUnitCapa)
			dt=z->data[m-nUnitCapa];
		else
			dt=0;
		data[m]=dt*d1_FOR;
	}
	/*计算移位后和原始序列重叠的部分*/
	for(m=nSizeOriginal-1;m>=nUnitCapa;m--)
	{
		dt=z->data[m];
		if (m>=nUnitCapa)
			dt2=z->data[m-nUnitCapa];
		else
			dt2=0;
		data[m]=dt*dFOR+dt2*d1_FOR;
	}
	/*计算原始序列且不和移位后重叠的部分*/
	for(m=0;m<nUnitCapa;m++)
	{
		if (m<nSizeOriginal)
			dt=z->data[m];
		else
			dt=0;
		data[m]=dFOR*dt;
	}
	//计算新的LOLP
	for(m=0;m<nLoad;m++)
	{
		if(m>=nSizeplus)break;
		if(nIsEENS==0)
			dLE+=data[m];
		else 
			dLE+=data[m]*(nLoad-m)*WPS_RELIABILITY_STEP;
	}

	free(data);
	return dLE;
}