#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include "blp.h"

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


blp* blp_create(int m,int n,int rn[],int cn[])
{
	blp* p;
	int cn111[]={1},i;
	/*allocate memory for blp object*/
	p=malloc(sizeof(blp));

	p->A=bmtx_create(m,n,rn,cn);
	p->a=bmtx_create(m,1,rn,cn111);
	p->b=bmtx_create(m,1,rn,cn111);
	p->c=bmtx_create(n,1,cn,cn111);
	p->d=bmtx_create(n,1,cn,cn111);
	p->e=bmtx_create(n,1,cn,cn111);
	p->x=bmtx_create(n,1,cn,cn111);
	p->xi=bmtx_create(n,1,cn,cn111);

	p->t=malloc(sizeof(int)*n);
	for(i=0;i<n;i++)p->t[i]=LPX_CV;
	
	/*allocate LPX problem object and initialize them*/
	p->lp=lpx_create_prob();
	lpx_add_rows(p->lp,p->A->rm);
	lpx_add_cols(p->lp,p->A->rn);
	return p;
	
}

void blp_free(blp* p)
{
	if(p==NULL)return;
	bmtx_free(p->A);
	bmtx_free(p->a);
	bmtx_free(p->b);
	bmtx_free(p->c);
	bmtx_free(p->d);
	bmtx_free(p->e);
	bmtx_free(p->x);
	bmtx_free(p->xi);
	free(p->t);
	lpx_delete_prob(p->lp);
	free(p);
}

/*set structural variable's kind*/
void blp_set_col_kind(blp* p,int m,int t)
{
	if(m<1||m>p->A->n)
		fault("blp_set_col_kind:index out of range");
	if(!(t==LPX_CV||t==LPX_IV01||t==LPX_IV))
		fault("blp_set_col_kind:unknow col type");
	p->t[m-1]=t;
}


void blp_load_c(blp* p)
{
	int i,j,bas,bn,nz=bmtx_get_nonzeros(p->c);
	mtx* m;
	for(i=1;i<=p->c->m;i++)
	{
		bn=bmtx_get_bk_cols(p->A,i);
		m=bmtx_get_bk_read(p->c,i,1);
		bas=bmtx_get_bk_row_base(p->c,i);
		if(mtx_is_valid(m))
		{
			for(j=1;j<=bn;j++)
			{
				lpx_set_obj_coef(p->lp,bas+j,mtx_get_element(m,j,1));
			}
		}
		else
		{
			for(j=1;j<=bn;j++)
			{
				lpx_set_obj_coef(p->lp,bas+j,0);
			}
		}
	}
}

int blp_load_get_bnd_type(int bl,int bu)
{
	int bt;
	if(bl!=0&&bu!=0)bt=LPX_DB;
	else if(bl!=0&&bu==0)bt=LPX_LO;
	else if(bl==0&&bu!=0)bt=LPX_UP;
	else if(bl==0&&bu==0)bt=LPX_FR;
	return bt;
}

int blp_equal_almost(double a,double b)
{
	return (fabs(a-b)<1e-6);
}

void blp_load_ab(blp* p)
{
	int i,j,bas,bl,bu,bt,bm;
	mtx *l,*u;
	double dl,du;
	for(i=1;i<=p->a->m;i++)
	{
		l=bmtx_get_bk_read(p->a,i,1);
		u=bmtx_get_bk_read(p->b,i,1);
		bl=mtx_is_valid(l);
		bu=mtx_is_valid(u);
		bas=bmtx_get_bk_row_base(p->a,i);
		bm=bmtx_get_bk_rows(p->A,i);
		for(j=1;j<=bm;j++)
		{
			bt=blp_load_get_bnd_type(bl,bu);
			dl=((bl!=0)?mtx_get_element(l,j,1):0);
			du=((bu!=0)?mtx_get_element(u,j,1):0);

			/*if lower boundary and upper boundary are equal almost,then set boundary-type to be LPX_FX*/
			if(bt==LPX_DB&&blp_equal_almost(dl,du)!=0)bt=LPX_FX;
			
			lpx_set_row_bnds(p->lp,bas+j,bt,dl,du);
			
			/*对于人工变量边界类型是LPX_FR的，打上标记，最后全部清除*/
/*
			if(bt==LPX_FR)
			{
				lpx_mark_row(p->lp,bas+j,1);
			}
*/
		}
	}
}

void blp_load_de_t(blp* p)
{
	int i,j,bas,bl,bu,bt,vt,bn,No_IV;
	mtx *l,*u;
	double dl,du;
	No_IV=0;
	for(i=1;i<=p->d->m;i++)
	{
		l=bmtx_get_bk_read(p->d,i,1);
		u=bmtx_get_bk_read(p->e,i,1);
		bl=mtx_is_valid(l);
		bu=mtx_is_valid(u);
		bas=bmtx_get_bk_row_base(p->d,i);
		vt=p->t[i-1];
		if(vt!=LPX_CV&&lpx_get_class(p->lp)!=LPX_MIP)
		{
			lpx_set_class(p->lp,LPX_MIP);
		}
		bn=bmtx_get_bk_cols(p->A,i);
		for(j=1;j<=bn;j++)
		{
			bt=blp_load_get_bnd_type(bl,bu);
			vt=p->t[i-1];
			dl=((bl!=0)?mtx_get_element(l,j,1):0);
			du=((bu!=0)?mtx_get_element(u,j,1):0);
			/*对于整数类型变量，将其边界整数化*/
			/*
			if(vt!=LPX_CV)
			{
				dl=ceil(dl);du=floor(du);
			}
			*/
			/*对于定义为0-1的变量，设定缺省边界，并且验证0-1*/
			if(vt==LPX_IV01)
			{
				if(bl==0)dl=0;
				if(bu==0)du=1;
				bt=LPX_DB;
				if(dl*(dl-1)!=0||du*(du-1)!=0)
				{
					fault("");
					/*error*/
				}
				/*LPX_IV01是自己添加的类型，所以要变成LPX认识的类型，即LPX_IV*/
				vt=LPX_IV;
				No_IV++;
			}
			if(bt==LPX_DB&&blp_equal_almost(dl,du)!=0)
			{
				bt=LPX_FX;
				if (vt==LPX_IV)
				{
					No_IV--;
				}
				/*对于固定值，强制指定为LPX_CV*/
				vt=LPX_CV;
			}
			if(lpx_get_class(p->lp)==LPX_MIP)lpx_set_col_kind(p->lp,bas+j,vt);
			lpx_set_col_bnds(p->lp,bas+j,bt,dl,du);
		}
	}
	if (No_IV<=0)lpx_set_class(p->lp,LPX_LP);
}


void blp_load_A(blp* p,ssm_treat_fun pstf)
{
	ssm* pssm=bmtx_create_ssm(p->A);
	//bmtx_clear_all(p->A);
	if(pstf!=NULL)
	{
		if(pssm==NULL)pssm=ssm_create(0);
		(*pstf)(pssm);
	}
	if(pssm!=NULL)
	{
		if(pssm->n>1)
			lpx_load_matrix(p->lp,pssm->n,pssm->rn,pssm->cn,pssm->a);
		ssm_free(pssm);
	}
}

int blp_get_col_base(blp* p,int n)
{
	return bmtx_get_bk_col_base(p->A,n);
}

int blp_get_row_base(blp* p,int m)
{
	return bmtx_get_bk_row_base(p->A,m);
}


//#define BLP_PRINT_KERNEL_TIME
int blp_solve(blp* p,ssm_treat_fun pstf)
{
	int ret,i,j,bas;
	mtx* x;
	double* px;
	int nMaxIterations;
	double dTimeout;
#ifdef BLP_PRINT_KERNEL_TIME
	double dt=0;
	clock_t t_start;
#endif
/*
	lpx_unmark_all(p->lp);
*/
	blp_load_A(p,pstf);
	blp_load_ab(p);
	/*对于约束类型为LPX_FR的，清除之，可以加快速度*/
/*
	lpx_del_items(p->lp);
*/
	lpx_get_class(p->lp);
	blp_load_de_t(p);
	blp_load_c(p);
	lpx_get_class(p->lp);

	/*todo:添加规划分析代码，清除上下限相等的变量，降低维度*/

	//printf("--------------------start--------------------\n");

	/*lpx_write_lpt(p->lp,"lpt.txt");*/

#ifdef BLP_PRINT_KERNEL_TIME
	t_start=clock();
#endif

	{
		nMaxIterations=lpx_get_int_parm(p->lp,LPX_K_ITLIM);
		dTimeout=lpx_get_real_parm(p->lp,LPX_K_TMLIM);
		lpx_set_int_parm(p->lp,LPX_K_ITLIM,-1);
		lpx_set_real_parm(p->lp,LPX_K_TMLIM,-1);
	}
	
	ret=lpx_simplex(p->lp);	
#ifdef BLP_PRINT_KERNEL_TIME
	dt+=(((double)(clock()-t_start))/CLOCKS_PER_SEC);
#endif
	if(lpx_get_status(p->lp)!=LPX_FEAS
		&&lpx_get_status(p->lp)!=LPX_OPT
		)
		return ret;

	/*保存非整数优化结果*/
	for(i=1;i<=p->x->m;i++)
	{
		bas=bmtx_get_bk_row_base(p->x,i);
		x=bmtx_get_bk_write(p->x,i,1);
		px=mtx_get_buf_write(x);
		for(j=1;j<=x->m;j++)
		{
			*px=lpx_get_col_prim(p->lp,bas+j);
			px++;
		}
	}
	p->obj=lpx_get_obj_val(p->lp);
	if(lpx_get_int_parm(p->lp,LPX_K_MSGLEV)>0)
		printf("\nobj=%12.4lf\n",p->obj);

	if(lpx_get_class(p->lp)==LPX_MIP)
	{
#ifdef BLP_PRINT_KERNEL_TIME
		t_start=clock();
#endif
		lpx_set_int_parm(p->lp,LPX_K_ITLIM,nMaxIterations);
		lpx_set_real_parm(p->lp,LPX_K_TMLIM,dTimeout);
		ret=lpx_integer(p->lp);

#ifdef BLP_PRINT_KERNEL_TIME
		dt+=(((double)(clock()-t_start))/CLOCKS_PER_SEC);
		printf("\t%7.3lf\t",dt);
#endif
		if(lpx_mip_status(p->lp)==LPX_I_OPT
			||lpx_mip_status(p->lp)==LPX_I_FEAS)
		{
			/*保存整数优化结果*/
			for(i=1;i<=p->xi->m;i++)
			{
				bas=bmtx_get_bk_row_base(p->xi,i);
				x=bmtx_get_bk_write(p->xi,i,1);
				px=mtx_get_buf_write(x);
				for(j=1;j<=x->m;j++)
				{
					*px=lpx_mip_col_val(p->lp,bas+j);
					px++;
				}
			}
			p->obji=lpx_mip_obj_val(p->lp);

			if(lpx_get_int_parm(p->lp,LPX_K_MSGLEV)>0)
				printf("\nobj=%12.4lf\n",p->obji);
		}
	}
	return ret;
}

/*Form GLPK format structures p->LP using the preset data in p*/
void blp_load_lp(blp* p,ssm_treat_fun pstf)
{
	blp_load_A(p,pstf);
	blp_load_ab(p);
	blp_load_de_t(p);
	blp_load_c(p);
}

/*load result from outside e.g. double *fs_result*/
void blp_load_result(blp* p,double *fs_result,double fs_fobjval)
{
	int i,j,bas,vt;
	mtx* x;
	double* px;

	/*保存非整数优化结果*/
	if(lpx_get_class(p->lp)==LPX_LP)
	{
		for(i=1;i<=p->x->m;i++)
		{
			bas=bmtx_get_bk_row_base(p->x,i);
			x=bmtx_get_bk_write(p->x,i,1);
			px=mtx_get_buf_write(x);
			for(j=1;j<=x->m;j++)
			{
				if (fs_result[bas+j]<TINY_0000001 && fs_result[bas+j]>-TINY_0000001)
					fs_result[bas+j]=0;
				*px=fs_result[bas+j];
				px++;
			}
		}
		p->obj=fs_fobjval;
	}

	if(lpx_get_class(p->lp)==LPX_MIP)
	{
		/*保存整数优化结果*/
		for(i=1;i<=p->xi->m;i++)
		{
			bas=bmtx_get_bk_row_base(p->xi,i);
			x=bmtx_get_bk_write(p->xi,i,1);
			px=mtx_get_buf_write(x);
			//Cplex有可能给出1.7e-17这样的数，所以对整数变量需要进行取整
			vt=p->t[i-1];
			if (vt==LPX_IV01||vt==LPX_IV)
			{
				for(j=1;j<=x->m;j++)
				{
					fs_result[bas+j]=floor(fs_result[bas+j]+0.499999);
				}
			}
			for(j=1;j<=x->m;j++)
			{
				if (fs_result[bas+j]<TINY_0000001 && fs_result[bas+j]>-TINY_0000001)
					fs_result[bas+j]=0;
				*px=fs_result[bas+j];
				px++;
			}
		}
		p->obji=fs_fobjval;
	}
}

///////////////////////以下函数将blp格式的优化模型改为能够用于Cplex的矩阵模型///////////////////////////
void blp_lp_malloc(blp* p,int *An,int *Am,int *nn,int **Ia,int **Ja,double **Ar,double **Qm_down,double **Qm_up,\
				   double **Xn_up,double **Xn_down,double **Cn,double **result,int **xstatus)
{
	int No_Am;
	int No_An;
	int No_nn;

	No_Am=lpx_get_num_rows(p->lp);
	No_An=lpx_get_num_cols(p->lp);
	No_nn=lpx_get_num_nz(p->lp);

	(*Am)=No_Am;
	(*An)=No_An;
	(*nn)=No_nn;

	(*Ia)=(int *)malloc(sizeof(int)*(No_nn+1));
	memset((*Ia),0,sizeof(int)*(No_nn+1));
	
	(*Ja)=(int *)malloc(sizeof(int)*(No_nn+1));
	memset((*Ja),0,sizeof(int)*(No_nn+1));
	
	(*Ar)=(double *)malloc(sizeof(double)*(No_nn+1));
	memset((*Ar),0,sizeof(double)*(No_nn+1));
	
	(*Qm_up)=(double *)malloc(sizeof(double)*(No_Am+1));
	memset((*Qm_up),0,sizeof(double)*(No_Am+1));
	
	(*Qm_down)=(double *)malloc(sizeof(double)*(No_Am+1));
	memset((*Qm_down),0,sizeof(double)*(No_Am+1));
	
	(*Xn_up)=(double *)malloc((No_An+1)*sizeof(double));
	memset((*Xn_up),0,sizeof(double)*(No_An+1));
	
	(*Xn_down)=(double *)malloc((No_An+1)*sizeof(double));
	memset((*Xn_down),0,sizeof(double)*(No_An+1));
	
	(*Cn)=(double *)malloc((No_An+1)*sizeof(double));
	memset((*Cn),0,sizeof(double)*(No_An+1));

	(*result)=(double *)malloc((No_An+1)*sizeof(double));
	memset((*result),0,sizeof(double)*(No_An+1));

	(*xstatus)=(int *)malloc(sizeof(int)*(No_An+1));
	memset((*xstatus),0,sizeof(int)*(No_An+1));
}

void blp_lp_get_data(blp* p,int An,int Am,int nn,int *Ia,int *Ja,double *Ar,double *Qm_down,double *Qm_up,\
				   double *Xn_up,double *Xn_down,double *Cn,int *xstatus)
{
	int i,j,No_element,bas,type;
	int * temp_ind;
	double * temp_val;

	temp_ind=(int *)malloc((An+1)*sizeof(int));
	memset((temp_ind),0,sizeof(int)*(An+1));
	temp_val=(double *)malloc((An+1)*sizeof(double));
	memset((temp_val),0,sizeof(double)*(An+1));
	
	//The matrix
	bas=0;
	for (i=1;i<=Am;i++)
	{
		No_element=lpx_get_mat_row(p->lp,i,temp_ind,temp_val);
		for (j=1;j<=No_element;j++)
		{
			Ia[j+bas]=i;
			Ja[j+bas]=temp_ind[j];
			Ar[j+bas]=temp_val[j];
		}
		bas=bas+No_element;
	}

	//obj_coefficient
	for (i=1;i<=An;i++)
	{
		Cn[i]=lpx_get_obj_coef(p->lp,i);
	}

	//The bound of rows
	for (i=1;i<=Am;i++)
	{
		type=lpx_get_row_type(p->lp,i);
		switch (type)
		{
			case LPX_DB:
				Qm_up[i]=lpx_get_row_ub(p->lp,i);
				Qm_down[i]=lpx_get_row_lb(p->lp,i);
				break;
			case LPX_LO:
				Qm_up[i]=INFINITE_1000000;
				Qm_down[i]=lpx_get_row_lb(p->lp,i);
				break;
			case LPX_UP:
				Qm_up[i]=lpx_get_row_ub(p->lp,i);
				Qm_down[i]=-INFINITE_1000000;
				break;
			case LPX_FX:
				Qm_up[i]=lpx_get_row_lb(p->lp,i);
				Qm_down[i]=lpx_get_row_lb(p->lp,i);
				break;
			default:
				Qm_up[i]=INFINITE_1000000;
				Qm_down[i]=-INFINITE_1000000;
				break;
		}
	}

	//The bound of rows
	for (i=1;i<=An;i++)
	{
		type=lpx_get_col_type(p->lp,i);
		switch (type)
		{
			case LPX_DB:
				Xn_up[i]=lpx_get_col_ub(p->lp,i);
				Xn_down[i]=lpx_get_col_lb(p->lp,i);
				break;
			case LPX_LO:
				Xn_up[i]=INFINITE_1000000;
				Xn_down[i]=lpx_get_col_lb(p->lp,i);
				break;
			case LPX_UP:
				Xn_up[i]=lpx_get_col_ub(p->lp,i);
				Xn_down[i]=-INFINITE_1000000;
				break;
			case LPX_FX:
				Xn_up[i]=lpx_get_col_lb(p->lp,i);
				Xn_down[i]=lpx_get_col_lb(p->lp,i);
				break;
			default:
				Xn_up[i]=INFINITE_1000000;
				Xn_down[i]=-INFINITE_1000000;
				break;
		}
	}

	//bound of variables LPX_IV or LPX_CV
	if (lpx_get_class(p->lp)==LPX_MIP)
	{
		for (i=1;i<=An;i++)
		{
			type=lpx_get_col_kind(p->lp,i);
			if(type==LPX_IV)
				xstatus[i]=1;
			else if(type==LPX_CV)
				xstatus[i]=0;
			else
			{fault("");}
		}
	}
	else
	{
		for (i=1;i<=An;i++)
		{
			xstatus[i]=0;
		}
	}
	free(temp_ind);
	free(temp_val);
}

void blp_lp_free(blp* p,int *Ia,int *Ja,double *Ar,double *Qm_down,double *Qm_up,\
				   double *Xn_up,double *Xn_down,double *Cn,double *result,int *xstaus)
{
	free(Ia);
	free(Ja);
	free(Ar);
	
	free(Qm_up);
	free(Qm_down);

	free(Xn_up);
	free(Xn_down);

	free(Cn);
	free(result);
	free(xstaus);
}

#ifdef _MEMINFO_ZHOUANSHI_
void lib_check_mem_leak()
{
LIBMEM* p;
LIBENV* env=lib_env_ptr();
printf("正在检测内存泄露...\n");
if(env->mem_count>0)
{
p=env->mem_ptr;
while(p->prev!=NULL)p=p->prev;
while(p!=NULL)
{
printf("	 %d\n",p->size);
p=p->next;
}
}
}
#endif

