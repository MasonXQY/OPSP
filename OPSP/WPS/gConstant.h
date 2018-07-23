  
#ifndef __G_CONSTANT_ZHANGNING_130307_H__
#define __G_CONSTANT_ZHANGNING_130307_H__

#define RT_WIND_SMLT 30/*������ģ�����*/
#define RT_WIND_FCST 31/*������ģ��Ԥ�����*/
#define RT_WIND_FERR  32/*������Ԥ����Ԥ��-ʵ�ʣ�*/

#define RT_WIND_STAT_PDF 33/*��糡���������ܶ�ͳ��*/
#define RT_WIND_STAT_CDF 34/*��糡�����ۻ�����ͳ��*/
#define RT_WIND_ALL_STAT_PDF 35/*����ܳ���ͳ��*/
#define RT_WIND_ALL_STAT_CDF 36/*����ܳ���ͳ��*/
#define RT_WIND_ALL_TYPICAL_SIMULATE 37/*���ģ����ͳ�������*/
#define RT_WIND_ALL_TYPICAL_THEORETICAL 38/*������۵��ͳ�������*/
#define RT_WIND_ALL_RESERVE_TYPICAL 39/*���Ԥ������Ⱥ�����������1��Ϊ�������������2��Ϊ���������������3��Ϊ�Ⱥɸ���������4��Ϊ�Ⱥɸ�����������EX4Ϊָʾ������0��ʾ�����������1��ʾƽ�������2��ʾ�����������*/
#define RT_WIND_ALL_RESERVE 40/*������ˮƽ�·�籸��������ǰ20���Ƿ��Ԥ��ֵ��ˮƽ���������ô�С��EX4Ϊָʾ������0��ʾ�����ã�1��ʾ������*/


/*��������״̬*/
#define US_IN_APP_ON 1/*ָ������*/
#define US_IN_APP_OFF 0/*ָ���ػ�*/

/*���򷵻�ֵ*/
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

#define WPS_DISTRIBUTION_INTERVAL 20 //���������ʷֲ��ĵ�����ĿǰĬ��Ϊ��Ϊ20�Ρ�
#define WPS_MIN_OUTPUT_CONFIDENCE 0.95 //��֤����ʱ�����Ŷ�
#define WPS_MAX_OUTPUT_CONFIDENCE 0.95 //��������ʱ�����Ŷ�
#define WPS_MIN_REGULATION_CONFIDENCE 0.99 //�����������������ʱ�����Ŷ�
#define WPS_MAX_REGULATION_CONFIDENCE 0.99 //�����������������ʱ�����Ŷ�
#define WPS_MIN_RESERVE_CONFIDENCE 0.99 //����󸺱���ʱ�����Ŷ�
#define WPS_MAX_RESERVE_CONFIDENCE 0.99 //����������õ����Ŷ�

#define WPS_POS_NEG_REGULATION 0.025 //�������������ж��ż���������ı�ϵͳ��Ȳ�ķ��ȳ�����һ��ֵ������Ϊ��������򷴵��塣

#define WPS_SIMULATION_TIMES 12/*��糡����ģ�������ÿ��ģ�������ô��Σ�ѡ�����е�����ӽ�ƽ������������������Ϊ����ģ�����ݣ�*/

#endif __G_CONSTANT_ZHANGNING_130307_H__