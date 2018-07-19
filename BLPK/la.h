#ifndef __MATRIX_H_ZHOUANSHI_005518__
#define __MATRIX_H_ZHOUANSHI_005518__
#if 1

#ifdef _DEBUG
#include "..\GLPK\include\glpk.h"
#define malloc umalloc
#define calloc ucalloc
//#define fopen ufopen
//#define fclose ufclose
#define free ufree
#endif
#endif
/*LINEAR ALGEBRA*/

/*test*/
void la_main(void);

double* la_zeros(double* a,int m,int n);
double* la_ones(double *a,double e,int m,int n);

double* la_unit(double *a,int n,double e);
double* la_rands(double *a,double k,int m,int n);


double* la_sub(double* b,double* a,int m,int n,
			   int m1,int m2,int n1,int n2,int* mb,int* nb);
double* la_diag(double* b ,double* a,int m,int n,int r,int *pmb,int *pnb);
double* la_triu_self(double *a,int m,int n,int r);
double* la_triu(double* b,double *a,int m,int n,int r);
double* la_tril_self(double *a,int m,int n,int r);
double* la_tril(double* b,double *a,int m,int n,int r);



void la_print(double *a,int m,int n,char *fmt);
void la_fprint(void *fp/*FILE* fp*/,double *a,int m,int n,char *fmt);





/*b=a*/
double* la_copy(double* b,const double* a,int m,int n);
/*c=a*b*/
double* la_multiply(double *c,double *a,double *b,int m,int k,int n);
/*c=a.b*/
double* la_dot_multiply(double *c,double *a,double *b,int m,int n);
/*a=a.b*/
double* la_dot_multiply_self(double *a,double *b,int m,int n);
/*d=a*b*c*/
double* la_multiply3(double *d,double *a,double *b,double *c,int m,int k,int l,int n);
/*e=a*b*c*d*/
double* la_multiply4(double *e,double *a,double *b,double *c,double *d,int m,int k,int l,int x,int n);
/*f=a*b*c*d*e*/
double* la_multiply5(double *f,double *a,double *b,double *c,double *d,double *e
					 ,int m,int k,int l,int x,int y,int n);
/*c=a+b*/
double* la_add(double* c,double* a,double* b,int m,int n);
/*a=a+b*/
double* la_add_self(double* a,double*b,int m,int n);
/*c=a-b*/
double* la_minus(double* c,double* a,double* b,int m,int n);
/*a=a-b*/
double* la_minus_self(double* a,double*b,int m,int n);
/*b=k*a*/
double* la_multiply_number(double* c,double k,double* a,int m,int n);
/*a=k*a*/
double* la_multiply_number_self(double k,double* a,int m,int n);
/*a=inv(a)*/
double* la_invert_self(double* a,int n);
/*b=inv(a)*/
double* la_invert(double* b,double* a,int n);
/*a=a'*/
double* la_transpose_self(double *a,int m,int n);
/*b=a'*/
double* la_transpose(double* b,double *a,int m,int n);


#endif //__MATRIX_H_ZHOUANSHI_005518__