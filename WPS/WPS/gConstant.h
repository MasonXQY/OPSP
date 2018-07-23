  
#ifndef __G_CONSTANT_ZHANGNING_130307_H__
#define __G_CONSTANT_ZHANGNING_130307_H__

#define RT_WIND_SMLT 30/*风电机组模拟出力*/
#define RT_WIND_FCST 31/*风电机组模拟预测出力*/
#define RT_WIND_FERR  32/*风电机组预测误差（预测-实际）*/

#define RT_WIND_STAT_PDF 33/*风电场出力概率密度统计*/
#define RT_WIND_STAT_CDF 34/*风电场出力累积概率统计*/
#define RT_WIND_ALL_STAT_PDF 35/*风电总出力统计*/
#define RT_WIND_ALL_STAT_CDF 36/*风电总出力统计*/
#define RT_WIND_ALL_TYPICAL_SIMULATE 37/*风电模拟典型出力曲线*/
#define RT_WIND_ALL_TYPICAL_THEORETICAL 38/*风电理论典型出力曲线*/
#define RT_WIND_ALL_RESERVE_TYPICAL 39/*风电预测误差峰谷荷正备用需求，1列为峰荷正备用需求，2列为峰荷正备用率需求，3列为谷荷负备用需求，4列为谷荷负备用率需求，EX4为指示变量，0表示正调峰情况，1表示平均情况，2表示反调峰情况，*/
#define RT_WIND_ALL_RESERVE 40/*各出力水平下风电备用率需求，前20列是风电预测值各水平下正负备用大小，EX4为指示变量，0表示正备用，1表示负备用*/


/*机组输入状态*/
#define US_IN_APP_ON 1/*指定开机*/
#define US_IN_APP_OFF 0/*指定关机*/

/*程序返回值*/
#define WPS_EXIT_SUCCESS 0
#define WPS_EXIT_FAILURE 1
#define WPS_EXIT_PARALACK 2
#define WPS_EXIT_PARAERROR 3
#define WPS_OPEN_OUTPUT_FILE_FAILURE 5
#define WPS_OPEN_INPUT_FILE_FAILURE 6

#define  WPS_RELIABILITY_STEP 1

#define TINY 0.000001
#define INFINITY 100000000

#define TRUE 1
#define FALSE 0

#define WPS_DISTRIBUTION_INTERVAL 20 //风电输出概率分布的点数，目前默认为分为20段。
#define WPS_MIN_OUTPUT_CONFIDENCE 0.95 //求保证出力时的置信度
#define WPS_MAX_OUTPUT_CONFIDENCE 0.95 //求最大出力时的置信度
#define WPS_MIN_REGULATION_CONFIDENCE 0.99 //求理论正反调峰出力时的置信度
#define WPS_MAX_REGULATION_CONFIDENCE 0.99 //求理论正反调峰出力时的置信度
#define WPS_MIN_RESERVE_CONFIDENCE 0.99 //求最大负备用时的置信度
#define WPS_MAX_RESERVE_CONFIDENCE 0.99 //求最大正备用的置信度

#define WPS_POS_NEG_REGULATION 0.025 //理论正反调峰判定门槛，如果风电改变系统峰谷差的幅度超过这一比值，则认为是正调峰或反调峰。

#define WPS_SIMULATION_TIMES 12/*风电场出力模拟次数（每次模拟均做这么多次，选择其中电量最接近平均电量的那组数据作为最终模拟数据）*/

#endif __G_CONSTANT_ZHANGNING_130307_H__