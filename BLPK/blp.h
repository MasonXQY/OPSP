#ifndef __BLP_H_ZHOUANSHI_005518__
#define __BLP_H_ZHOUANSHI_005518__
/*
�޸�GLPK���룺
��ӷ����ڴ��ĳ�����붨λ��Ϣ


-----------------------------------------
�޸�LIBMEM�ṹ�壬������´���
#ifdef _DEBUG
#define _MEMINFO_ZHOUANSHI_
#endif
#ifdef _MEMINFO_ZHOUANSHI_
	  char info[65];
#endif

-----------------------------------------
�޸�umalloc��������
#ifdef _MEMINFO_ZHOUANSHI_
#define glp_lib_umalloc umalloc
#undef umalloc
#define umalloc(x) _umalloc(x,__FILE__, __LINE__)

void *_umalloc(int size,char* sf,int nl);
#else
void *umalloc(int size);
#endif

-----------------------------------------
�޸�umalloc������ʵ��
#ifdef _MEMINFO_ZHOUANSHI_
void *_umalloc(int size,char* sf,int nl)
#else
void *umalloc(int size)
#endif
������������

��Ӵ���
#ifdef _MEMINFO_ZHOUANSHI_
	  _snprintf(desc->info,64,"%s(%d)[%d]",sf,nl,size-size_of_desc);
#endif

-----------------------------------------
���Ӻ������ú����Ὣ���з������ڴ棬��û���ͷŵĴ���λ���оٳ���
#ifdef _MEMINFO_ZHOUANSHI_
void lib_check_mem_leak()
{
	LIBMEM* p;
	LIBENV* env=lib_env_ptr();
	if(env->mem_count>0)
	{
		p=env->mem_ptr;
		while(p->prev!=NULL)p=p->prev;
		while(p!=NULL)
		{
			printf("%s\n",p->info);
			p=p->next;
		}
	}
}
#else
#define lib_check_mem_leak
#endif
-----------------------------------------
�޸�spx_simplex�������ڶ�ż�����ͳ�����ֵ���ȶ������ʱ������spx->tm_lim=60;
�޸�mip_driver������
		 case LPX_E_TMLIM:
			if (tree->msg_lev >= 1)
               print("numerical instability (primal simplex, phase I OR II)"
			   ", TIMEOUT FORCED");
			lpx_set_real_parm(tree->lp,LPX_K_TMLIM,-1);
			goto fath;
            break;
-----------------------------------------
�޸�double lib_get_time(void)
{
	return ((double)clock())/CLOCKS_PER_SEC;
}
*/

#include "..\GLPK\include\glpk.h"
#include "mtx.h"
#define LPX_IV01          162   /* 0-1 integer variable */
#define INFINITE_1000000 1000000  /*Very large number*/
#define TINY_0000001 0.000001  /*Very small number*/


/*
�ֿ����
                ��  ��
                ��x1��
  min [c1,c2,c3]��x2��
                ��x3��
                ��  ��
  s.t.
    ��  ��  ��           ����  ��  ��  ��
    ��a1��  ��A11 A12 A12����x1��  ��b1��
    ��a2���ܩ�A11 A12 A12����x2���ܩ�b2��
    ��  ��  ��           ����x3��  ��  ��
    ��  ��  ��           ����  ��  ��  ��
    ��  ��  ��  ��  ��  ��
    ��d1��  ��x1��  ��e1��
    ��d2���ܩ�x2���ܩ�e2��
    ��d3��  ��x3��  ��e3��
    ��  ��  ��  ��  ��  ��

���������
 	LPX* lp;
 	�����ֿ����col_blocks��3
 	ÿ���ֿ���������ͣ�int[3] cols_type={t1,t2,t3};//LPX_CV:ȱʡ;LPX_IV01:0-1;LPX_IV:����
 	ÿ���ֿ������ά�ȣ�int[3] cols_size={sizeof(x1),sizeof(x2),sizeof(x3)}
 	ÿ���ֿ������Ŀ�����棺double*[3] coef={&c1,&c2,&c3}//null��ʾ0
 	ÿ���ֿ���������ޣ�double*[3] col_lower_bound={&d1,&d2,&d3}//null��ʾ������
 	ÿ���ֿ���������ޣ�double*[3] col_upper_bound={&e1,&e2,&e3}//null��ʾ������;
						���޺�����ָ����ͬ����Ϊ��ʽԼ��
 	Լ�������ֿ����row_blocks��2
 	ÿ���ֿ�Լ����ά�ȣ�int[2] rows_size={sizeof(a1),sizeof(a2)}
 	ÿ���ֿ�Լ�����������ޣ�double*[2] row_lower_bound={&a1,&a2 }//null��ʾ������
 	ÿ���ֿ�Լ�����������ޣ�double*[2] row_upper_bound={&b1,&b2 }//null��ʾ������;
							���޺�����ָ����ͬ����Ϊ��ʽԼ��
 	A����double*[2][3] A={{&A11,&A12,&A13},{&A21,&A22,&A23}}//null��ʾ0

*/

typedef struct _blp
{
	LPX* lp;/*LPX pointer*/
	bmtx *A/*constraint matrix*/
		,*a/*lower boundary of auxiliary variable*/
		,*b/*upper boundary of auxiliary variable*/
		,*c/*coefficient of the objective function*/
		,*d/*lower boundary of structural variable*/
		,*e/*upper boundary of structural variable*/
		,*x/*�������Ż����*/
		,*xi/*�����Ż����*/
		;
	int *t/* x(i)�����ͣ�LPX_CV,LPX_IV,LPX_IV01 */;
	double obj/*�������Ż�Ŀ��ֵ*/
		,obji/*�����Ż�Ŀ��ֵ*/
		;
}blp;

/*blp test function entrance*/
void blp_main(void);

/*create blocked linear or mixed integer linear programming object*/
blp* blp_create(int m,int n,int rn[],int cn[]);

/*destory the programming object*/
void blp_free(blp* p);

/*set the kind of structural variables*/
void blp_set_col_kind(blp* p,int m,int t);

/*slove the programming*/
int blp_solve(blp* p,ssm_treat_fun pstf);

/*Form GLPK format structures p->LP using the preset data in p*/
void blp_load_lp(blp* p,ssm_treat_fun pstf);

/*load result from outside e.g. double *fs_result*/
void blp_load_result(blp* p,double *fs_result,double fs_fobjval);


void blp_lp_malloc(blp* p,int *An,int *Am,int *nn,int **Ia,int **Ja,double **Ar,double **Qm_down,double **Qm_up,\
				   double **Xn_up,double **Xn_down,double **Cn,double **result,int **xstatus);

void blp_lp_get_data(blp* p,int An,int Am,int nn,int *Ia,int *Ja,double *Ar,double *Qm_down,double *Qm_up,\
				   double *Xn_up,double *Xn_down,double *Cn,int *xstatus);

void blp_lp_free(blp* p,int *Ia,int *Ja,double *Ar,double *Qm_down,double *Qm_up,\
				   double *Xn_up,double *Xn_down,double *Cn,double *result,int *xstatus);

/**/
int blp_get_col_base(blp* p,int n);
int blp_get_row_base(blp* p,int m);

#ifdef _DEBUG
#define _MEMINFO_ZHOUANSHI_
#endif

#ifdef _MEMINFO_ZHOUANSHI_
/*check memory leak,it will print all memory block those had not freeed*/
void lib_check_mem_leak();
#else
#define lib_check_mem_leak()
#endif
#endif //__BLP_H_ZHOUANSHI_005518__