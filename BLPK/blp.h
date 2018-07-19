#ifndef __BLP_H_ZHOUANSHI_005518__
#define __BLP_H_ZHOUANSHI_005518__
/*
修改GLPK代码：
添加分配内存块的程序代码定位信息


-----------------------------------------
修改LIBMEM结构体，添加如下代码
#ifdef _DEBUG
#define _MEMINFO_ZHOUANSHI_
#endif
#ifdef _MEMINFO_ZHOUANSHI_
	  char info[65];
#endif

-----------------------------------------
修改umalloc函数定义
#ifdef _MEMINFO_ZHOUANSHI_
#define glp_lib_umalloc umalloc
#undef umalloc
#define umalloc(x) _umalloc(x,__FILE__, __LINE__)

void *_umalloc(int size,char* sf,int nl);
#else
void *umalloc(int size);
#endif

-----------------------------------------
修改umalloc函数的实现
#ifdef _MEMINFO_ZHOUANSHI_
void *_umalloc(int size,char* sf,int nl)
#else
void *umalloc(int size)
#endif
。。。。。。

添加代码
#ifdef _MEMINFO_ZHOUANSHI_
	  _snprintf(desc->info,64,"%s(%d)[%d]",sf,nl,size-size_of_desc);
#endif

-----------------------------------------
增加函数，该函数会将所有分配了内存，但没有释放的代码位置列举出来
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
修改spx_simplex函数，在对偶单纯型出现数值不稳定的情况时，增加spx->tm_lim=60;
修改mip_driver函数，
		 case LPX_E_TMLIM:
			if (tree->msg_lev >= 1)
               print("numerical instability (primal simplex, phase I OR II)"
			   ", TIMEOUT FORCED");
			lpx_set_real_parm(tree->lp,LPX_K_TMLIM,-1);
			goto fath;
            break;
-----------------------------------------
修改double lib_get_time(void)
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
分块矩阵
                ┌  ┐
                │x1│
  min [c1,c2,c3]│x2│
                │x3│
                └  ┘
  s.t.
    ┌  ┐  ┌           ┐┌  ┐  ┌  ┐
    │a1│  │A11 A12 A12││x1│  │b1│
    │a2│≤│A11 A12 A12││x2│≤│b2│
    │  │  │           ││x3│  │  │
    └  ┘  └           ┘└  ┘  └  ┘
    ┌  ┐  ┌  ┐  ┌  ┐
    │d1│  │x1│  │e1│
    │d2│≤│x2│≤│e2│
    │d3│  │x3│  │e3│
    └  ┘  └  ┘  └  ┘

输入参数：
 	LPX* lp;
 	变量分块个数col_blocks：3
 	每个分块变量的类型：int[3] cols_type={t1,t2,t3};//LPX_CV:缺省;LPX_IV01:0-1;LPX_IV:整数
 	每个分块变量的维度：int[3] cols_size={sizeof(x1),sizeof(x2),sizeof(x3)}
 	每个分块变量的目标增益：double*[3] coef={&c1,&c2,&c3}//null表示0
 	每个分块变量的下限：double*[3] col_lower_bound={&d1,&d2,&d3}//null表示无下限
 	每个分块变量的上限：double*[3] col_upper_bound={&e1,&e2,&e3}//null表示无上限;
						下限和上限指针相同，则为等式约束
 	约束条件分块个数row_blocks：2
 	每个分块约束的维度：int[2] rows_size={sizeof(a1),sizeof(a2)}
 	每个分块约束条件的下限：double*[2] row_lower_bound={&a1,&a2 }//null表示无下限
 	每个分块约束条件的上限：double*[2] row_upper_bound={&b1,&b2 }//null表示无上限;
							下限和上限指针相同，则为等式约束
 	A矩阵：double*[2][3] A={{&A11,&A12,&A13},{&A21,&A22,&A23}}//null表示0

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
		,*x/*非整数优化结果*/
		,*xi/*整数优化结果*/
		;
	int *t/* x(i)的类型，LPX_CV,LPX_IV,LPX_IV01 */;
	double obj/*非整数优化目标值*/
		,obji/*整数优化目标值*/
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