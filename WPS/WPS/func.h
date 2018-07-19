#include "gConstant.h"

//WPOS_sub/////////////////////////////////////////////////////////////////////////////
void app_err(int nErrorCode,char *fmt, ...);

void app_stat_wind_typical_curve(app* me,double *daLoad,double *dapower,int nCount,\
	double *dPosRegulation,double *dNegRegulation,double * dFlatRegulation,\
	double *dPosRegulationPro,double *dNegRegulationPro,double *dFlatRegulationPro);

double app_stat_load_peak_valley(double *daDailyLoad,int nCount,int nType);

void app_stat_wind_daily_curve(app* me, double *dAverageLoad,double *dapower,int nCount,double dPowerMaxall,\
	double *dPosRegulation,double *dNegRegulation,double *dMaxOutput,\
	double *dMinOutput,double *dAverageOutput);

void	app_stat_wind_variation(double *dapower,double *davariation,double *daLoad,int nCount, double dPowerMaxall,\
	double *daPos,double *daNeg,double *daParameter);

void	app_stat_wind_reserve(double *davariation,int nCount, double* dMaxPos,double *dMaxNeg,double *dmeanPos, double *dmeanNeg);

void app_out_wind_power(app* me);

void app_out_wind_power_stat(app* me);

void app_out_wind_fpower(app* me);

void app_out_wind_ferror(app* me);

void app_out_wind_fpower_stat( app* me);

void app_out_wind_power_stat_month( app* me, mtx* pmWindPower[WPOS_SIMULATION_TIMES],int nStartDate, int nEndDate);

void app_out_wind_fpower_stat_month( app* me, mtx* pmWindPower[WPOS_SIMULATION_TIMES],mtx* pmWindfPower[WPOS_SIMULATION_TIMES],int nStartDate, int nEndDate);

//pa//////////////////////////////////////////////////////////////////////////////
pa* pa_create();

void pa_clear(pa *p,int bFreeItem);

void pa_free(pa *p,int bFreeItem);

void pa_reallocate(pa* p,int cap);

void* pa_get(pa *p,int idx/*base 0*/);

void pa_set(pa *p,int idx/*base 0*/,void* item);

void pa_add(pa *p,void* item);

//app_cur//////////////////////////////////////////////////////////////////////////////
unit_para* app_cur_unit_para(unit* pUnit,int nDate);
unit_power* app_cur_unit_power(unit* pUnit,int nDate);
unit_power* app_cur_unit_fpower(unit* pUnit,int nDate);
wind_para* app_cur_wind_para(unit* pUnit,int nDate);
wind_zone_para* app_cur_windzone_para(wind_zone* pWindZone,int nDate);
load* app_cur_load(app* me,int nDate);
regulation* app_cur_regulationpro(app* me,int nDate);
//iDate//////////////////////////////////////////////////////////////////////////////
int idate_is_leap(int nYear);

int idate_days_of_year(int iYear);

int idate_days_of_month(int iYear,int iMonth);

int idate_day_idx_of_year(int iDate);

int idate_week_idx_of_year(int iDate);

int idate_month_idx_of_year(int iDate);

int idate_absolute_days(int nDate);

int idate_days_between(int iDate1,int iDate2);

int idate_next_month(int iDate);

int idate_next_year(int iDate);

int idate_first_day_of_month(int iDate);

int idate_next_day(int iDate);

int idate_is_valid(int iDate);

int idate_identify_double(double d);

int idate_first_day_of_year(int iDate);

//main//////////////////////////////////////////////////////////////////////////////

int util_double2int(double d);

//Wind//////////////////////////////////////////////////////////////////////////////

double imaths_gen_normal(double u,double g, int *refresh);

void app_generate_wind_power(app *me);

void Generate_wind_power_month(int nDateStart,int nDateEnd,app *me,mtx* pmWindPower, int *nRefresh);

void Generate_wind_forecasterror_month(int nDateStart,int nDateEnd,app *me,mtx* pmWindPower,mtx* pmWindFPower, int *nRefresh);

double Generate_wind_power_Nu(double x, double c, double k, double theta);

double Generate_wind_forecasterror_Nu(double theta);

double Generate_wind_power_curve(double dWindSP, double dCutinSP, double dRatedSP, double dCutoutSP);

double Generate_wind_power_curve_cubic(double dWindSP, double dCutinSP, double dRatedSP, double dCutoutSP);

double Generate_wind_power_curve_linear(double dWindSP, double dCutinSP, double dRatedSP, double dCutoutSP);

void Generate_wind_power_SetWindPower(int nDateStart,int nDateEnd,app *me,mtx* pmWindPower);

void Generate_wind_power_SetWindFPower(int nDateStart,int nDateEnd,app *me,mtx* pmWindFPower);

void Generate_wind_power_StatisticsWindPower(app *me);

double calc_WindCC_on_LOLPorEENS(mtx* z,double dLOLPorEENS, double dFor,int nLoad, double dmax, int nIsEENS);

double calc_LOLPorEENS(mtx* z,double dCapa, double dFor,int nLoad,int nIsEENS);

