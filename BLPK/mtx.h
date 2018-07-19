#ifndef __MATRIX_MTX_ZHOUANSHI_005518_H__
#define __MATRIX_MTX_ZHOUANSHI_005518_H__
#include "la.h"
#ifdef _DEBUG
#define _SAFE_CHECK_
#endif
typedef struct _ssm/*simple sparse matrix*/
{
	int n,*rn,*cn;
	double *a;
	int cap;/*分配的内存容量*/
}ssm;


typedef struct _mtx/*normal matrix*/
{
	int m,n;
	int lock;/*lock dimensions*/
	/*初始状态并不申请数据内存空间*/
	double* data;
	int buf_len;

#define MTX_NORMAL	0
#define MTX_SSM		1
	int t;
	ssm* pssm;
}mtx;

typedef struct _bmtx
{
	/*一旦创建，所有的维度信息都不能更改*/
	int m,n;/*rows and cols of block matrix*/
	int *bm,*bn;
	int *sbm,*sbn;
	int rm,rn;/*real rows and cols of real matrix*/

	/*初始状态并不创建mtx对象*/
	mtx** bks;/*mtx blocks*/
}bmtx;




/*------------------------------*------------------------------*/

/*create simple sparse matrix object*/
ssm* ssm_create(int n);
/*free the simple sparse matrix object*/
void ssm_free(ssm* p);
/*add element*/
void ssm_add(ssm* p,int r,int c,double a);
typedef void (*ssm_treat_fun)(ssm* p);
/*初始化矩阵为k倍的单位阵*/
ssm* ssm_create_unit(int n,double k);
ssm* ssm_create_diag(int n,double* data);
ssm* ssm_create_from_mtx(mtx* a);
/*------------------------------*------------------------------*/

mtx* mtx_t(int i,double* d,int m,int n);

double mtx_sum_element(mtx* p);
/*check the matrix p is valid or not
if data buffer had not been allocated or buffer size less then matrix size*/
int mtx_is_valid(mtx* p);
int mtx_is_invalid(mtx* p);
/*clear the matrix p's data buffer*/
void mtx_clear_buf(mtx* p);
void bmtx_clear_all(bmtx* p);

/*free the matrix object*/
void mtx_free(mtx* p);
/*get a safe data buffer for write*/
double* mtx_get_buf_write(mtx* p);
/*get a safe data buffer for read*/
double* mtx_get_buf_read(mtx* p);
/*set new dimensions of the matrix p*/
mtx* mtx_resize(mtx* r,int m,int n);
/*create new matrix object*/
mtx* mtx_create(int m,int n,int lock);
/*print the matrix to stdout or FILE*/
void mtx_print(mtx* a);
void mtx_fprint(void* fp,mtx* a);
/*r=a*/
mtx* mtx_copy(mtx* r,mtx* a);
mtx* mtx_cover(mtx* r,mtx* a);
mtx* mtx_copy_data(mtx* r,double* data,int m,int n);

/*r=diag(a,r)*/
mtx* mtx_diag(mtx* r,mtx* a,int k);
/*r=a+b*/
mtx* mtx_add(mtx* r,mtx* a,mtx* b);
/*a=a+b*/
mtx* mtx_add_self(mtx* a,mtx* b);
/*r=a-b*/
mtx* mtx_minus(mtx* r,mtx* a,mtx* b);
/*a=a-b*/
mtx* mtx_minus_self(mtx* a,mtx* b);
/*r=a.*b*/
mtx* mtx_dot_mul(mtx* r,mtx* a,mtx* b);
/*a=a.*b*/
mtx* mtx_dot_mul_self(mtx* a,mtx* b);
/*r=a*b*/
mtx* mtx_mul(mtx* r,mtx* a,mtx* b);
mtx* mtx_mul_some(mtx* r,...);

/*r=a*b*c*/
mtx* mtx_mul3(mtx* r,mtx* a,mtx* b,mtx* c);
/*r=a*b*c*d*/
mtx* mtx_mul4(mtx* r,mtx* a,mtx* b,mtx* c,mtx* d);
/*r=a*b*c*d*e*/
mtx* mtx_mul5(mtx* r,mtx* a,mtx* b,mtx* c,mtx* d,mtx* e);
/*r=k*a*/
mtx* mtx_mul_number(mtx* r,double k,mtx* a);
/*a=k*a*/
mtx* mtx_mul_number_self(double k,mtx* a);
/*r=inv(a)*/
mtx* mtx_invert(mtx* r,mtx* a);
mtx* mtx_invert_self(mtx* a);
/*r=a'*/
mtx* mtx_transpose(mtx* r,mtx* a);
/*a=a'*/
mtx* mtx_transpose_self(mtx* a);
mtx* mtx_sub(mtx* r,mtx* a,int m1,int m2,int n1,int n2);
mtx* mtx_tril(mtx* r,mtx* a,int k);
mtx* mtx_tril_self(mtx* a,int k);
mtx* mtx_triu(mtx* r,mtx* a,int k);
mtx* mtx_triu_self(mtx* a,int k);
/*初始化所有元素为0*/
mtx* mtx_ini_zeros(mtx* r,int m,int n);
/*初始化矩阵为k倍的单位阵*/
mtx* mtx_ini_unit(mtx* r,int m,double k);
/*初始化所有元素为k*/
mtx* mtx_ini_ones(mtx* r,int m,int n,double k);
/*初始化所有元素为0-k的随机数*/
mtx* mtx_ini_rands(mtx* r,int m,int n,double k);
mtx* mtx_ini_diag(mtx* r,int m,double* data);

double mtx_get_element(mtx* a,int m,int n);
void mtx_set_element(mtx* a,int m,int n,double e);

int mtx_get_nonzeros(mtx* a);


mtx* mtx_ORPlusConvolution(mtx* z,mtx* x,mtx* y);
mtx* mtx_ORPlusConvolution_self(mtx* z,mtx* x);
/*------------------------------*------------------------------*/


/*create block matrix object*/
bmtx* bmtx_create(int m,int n,int rows_bks[],int cols_bks[]);
/*get block*/
mtx* bmtx_get_bk_read(bmtx* p,int m,int n);
mtx* bmtx_get_bk_write(bmtx* p,int m,int n);
void bmtx_clear_bk(bmtx* p,int m,int n);
//double* bmtx_get_bk_buf_write(bmtx* p,int m,int n);
//double* bmtx_get_bk_buf_read(bmtx* p,int m,int n);
int bmtx_get_bk_rows(bmtx* p,int m);
int bmtx_get_bk_row_base(bmtx* p,int m);
int bmtx_get_bk_cols(bmtx* p,int n);
int bmtx_get_bk_col_base(bmtx* p,int n);
int bmtx_get_nonzeros(bmtx* p);
ssm* bmtx_create_ssm(bmtx* p);
void bmtx_free(bmtx* p);
void bmtx_set_bk_ssm(bmtx* p,int m,int n,ssm* pssm);

#endif //__MATRIX_MTX_ZHOUANSHI_005518_H__