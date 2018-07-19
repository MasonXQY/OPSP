#include "stdafx.h"


/*�ж��Ƿ�����*/
int idate_is_leap(int nYear)
{
	/*if(nYear%400==0||(nYear%4==0&&nYear%100!=0))�ж��ǲ�������
		return 1;
	else*/
		return 0;	
}
/*���ָ�����µ�����*/
int idate_days_of_year(int iYear)
{
	return 365+idate_is_leap(iYear);
}
/*���ָ�����µ�����*/
int idate_days_of_month(int iYear,int iMonth)
{
	switch(iMonth) 
	{
	case 1:case 3:case 5:case 7:case 8:case 10:case 12:return 31;
	case 4:case 6:case 9:case 11:return 30;
	case 2:return 28+idate_is_leap(iYear);
	default:return -1;
	}
}
/*��õ�ǰ�����������е�������1��1��Ϊ��һ��*/
int idate_day_idx_of_year(int iDate)
{
	int day=iDate%100,month=(iDate/100)%100,year=iDate/10000,sum=0;
	switch(month)/*�ȼ���ĳ����ǰ�·ݵ�������*/
	{
	case 1:sum=0;break;
	case 2:sum=31;break;
	case 3:sum=59;break;
	case 4:sum=90;break;
	case 5:sum=120;break;
	case 6:sum=151;break;
	case 7:sum=181;break;
	case 8:sum=212;break;
	case 9:sum=243;break;
	case 10:sum=273;break;
	case 11:sum=304;break;
	case 12:sum=334;break;
	}
	sum=sum+day;/*�ټ���ĳ�������*/
	if(month>2&&idate_is_leap(year)!=0)/*�ж��ǲ�������*/
		sum+=1;
	return sum;
}
/*��ø����ڵ���һ�����ǵڼ�������*/
int idate_week_idx_of_year(int iDate)
{
	return (idate_day_idx_of_year(iDate)-1)/7+1;
}
/*��ø����ڵ��·�*/
int idate_month_idx_of_year(int iDate)
{
	return (iDate/100)%100;
}

/*��õ�ǰ����20000101��������20000101Ϊ��һ��*/
int idate_absolute_days(int nDate)
{
	int sum=0;
	int i;
	for(i=2000;i<nDate/10000;i++)sum+=idate_days_of_year(i);
	sum+=idate_day_idx_of_year(nDate);
	return sum;
}

/*�����������֮��ļ����ͬһ����Ϊ0*/
int idate_days_between(int iDate1,int iDate2)
{
	return idate_absolute_days(iDate2)-idate_absolute_days(iDate1);
}

/*��ø��·ݵ���һ���·ݵĵ�һ��*/
int idate_next_month(int iDate)
{
	int iYearMonth=iDate/100;
	return (iYearMonth/100+(iYearMonth%100)/12)*10000+((iYearMonth%100)%12+1)*100+1;
}

/*���iDate�ĵڶ���ĵ�һ��*/
int idate_next_year(int iDate)
{
	int iYear=iDate/10000;
	return iYear*10000+10101;
}

/*��ø��·ݵĵ�һ��*/
int idate_first_day_of_month(int iDate)
{
	int iYearMonth=iDate/100;
	return iYearMonth*100+1;
}

/*��ø���ĵ�һ��*/
int idate_first_day_of_year(int iDate)
{
	return (iDate/10000)*10000+101;
}

/*��ø����ڵ���һ��*/
int idate_next_day(int iDate)
{
	int iYear=iDate/10000;
	int	iMonth=idate_month_idx_of_year(iDate);
	int iday=iDate%100;
	if (iday==idate_days_of_month(iYear,iMonth))
		return idate_next_month(iDate);
	else
		return iDate+1;
}

/*�ж������Ƿ���Ч*/
int idate_is_valid(int iDate)
{
	int i=iDate/10000;
	if(i>3000||i<1000)return 0;
	i=(iDate/100)%100;
	if(i>12||i<1)return 0;
	i=iDate%100;
	if(i<1||i>idate_days_of_month(iDate/10000,(iDate/100)%100))return 0;
	return 1;
}

int idate_identify_double(double d)
{
	int i;
	if(d<1000)return 0;
	else if(d<100000)d=d*10000+0.5;
	else if(d<10000000)d=d*100+0.5;
	else if(d<1000000000)d=d  +0.5;
	else return 0;
	i=(int)d;
	if((i%100)==0)i+=1;
	if(((i/100)%100)==0)i+=100;
	if(idate_is_valid(i)==0)return 0;
	else return i;
}