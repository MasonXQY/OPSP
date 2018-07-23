#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <LIMITS.H>
#include "la.h"
#include "mtx.h"

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


/*------------------------------*------------------------------*/

/*create simple sparse matrix object*/
ssm* ssm_create(int n)
{
	ssm* p;
#ifdef _SAFE_CHECK_
	if(n<0)
		fault("ssm_create:bad matrix dimension,n=%d",n);
#endif
	n++;
	p=malloc(sizeof(ssm));
	p->n=0;
	if(n/5>1000)n=(int)(n*1.2);
	else n+=1000;
	p->cn=malloc(sizeof(int)*n);
	p->rn=malloc(sizeof(int)*n);
	p->a=malloc(sizeof(double)*n);
	p->cap=n;
	*p->cn=5518;
	*p->rn=5518;
	*p->a=5518;
	return p;
}
/*free the simple sparse matrix object*/
void ssm_free(ssm* p)
{
	if(p==NULL)return;
	free(p->a);
	free(p->rn);
	free(p->cn);
	free(p);
}

void ssm_add(ssm* p,int r,int c,double a)
{
	void* pt;
	ssm* ps;
	p->n++;
	if(p->n>=p->cap)
	{
		ps=ssm_create(p->n);
		memcpy(ps->a,p->a,sizeof(double)*p->n);
		memcpy(ps->rn,p->rn,sizeof(int)*p->n);
		memcpy(ps->cn,p->cn,sizeof(int)*p->n);
		pt=p->a;p->a=ps->a;ps->a=pt;
		pt=p->rn;p->rn=ps->rn;ps->rn=pt;
		pt=p->cn;p->cn=ps->cn;ps->cn=pt;
		p->cap=ps->cap;
		ssm_free(ps);
	}
	p->cn[p->n]=c;
	p->rn[p->n]=r;
	p->a[p->n]=a;
}

/*初始化矩阵为k倍的单位阵*/
ssm* ssm_create_unit(int n,double k)
{
	int i;
	ssm* p=ssm_create(n);
	for(i=1;i<=n;i++)ssm_add(p,i,i,k);
	return p;
}

ssm* ssm_create_diag(int n,double* data)
{
	int i;
	ssm* p=ssm_create(n);
	for(i=1;i<=n;i++)ssm_add(p,i,i,data[i-1]);
	return p;
}
ssm* ssm_create_from_mtx(mtx* a)
{
	int i,j;
	int nz=mtx_get_nonzeros(a);
	ssm* p=ssm_create(nz);
	for(i=1;i<=a->m;i++)
		for(j=1;j<=a->n;j++)
			ssm_add(p,i,j,mtx_get_element(a,i,j));
	return p;
}

/*------------------------------*------------------------------*/

static mtx __mtx_[]={
	{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
	,{0,0,1,NULL,INT_MAX}
};

mtx* mtx_t(int i,double* d,int m,int n)
{
	mtx* p=__mtx_+i;
#ifdef _SAFE_CHECK_
	if(d==NULL||i<0||i>9)
		fault("");
#endif
	p->m=m;
	p->n=n;
	p->data=d;
	return p;
}

/*check the matrix p is valid or not
if data buffer had not been allocated or buffer size less then matrix size*/
int mtx_is_valid(mtx* p)
{
	if(p==NULL)return 0;
	if((p->data==NULL&&p->buf_len!=0)||p->buf_len<p->m*p->n)return 0;/*false*/
	return 1;/*true*/
}
int mtx_is_invalid(mtx* p)
{
	return mtx_is_valid(p)!=0?0:1;
}

/*clear the matrix p's data buffer*/
void mtx_clear_buf(mtx* p)
{
	if(p==NULL)return;
	p->buf_len=0;
	if(p->data!=NULL)
	{
		free(p->data);
		p->data=NULL;
	}
	if(p->pssm!=NULL)
	{
		ssm_free(p->pssm);
		p->pssm=NULL;
	}
}
/*free the matrix object*/
void mtx_free(mtx* p)
{
	if(p==NULL)return;
	mtx_clear_buf(p);
	free(p);
}

/*get a safe data buffer for write*/
double* mtx_get_buf_write(mtx* p)
{
#ifdef _SAFE_CHECK_

	if(p->t!=MTX_NORMAL)
		fault("unsupported!!");
#endif
	if(mtx_is_invalid(p))
	{
		mtx_clear_buf(p);
		p->buf_len=p->m*p->n;
		if(p->buf_len>0)
			p->data=malloc(sizeof(double)*p->buf_len);
	}
	return p->data;
}

/*get a safe data buffer for read*/
double* mtx_get_buf_read(mtx* p)
{
#ifdef _SAFE_CHECK_
	if(p->t!=MTX_NORMAL)
		fault("unsupported!!");
	if(mtx_is_invalid(p))
		fault("mtx_get_buf_read:data buffer invalid");
#endif
	return p->data;
}

/*set new dimensions of the matrix p*/
mtx* mtx_resize(mtx* r,int m,int n)
{
#ifdef _SAFE_CHECK_
	if(m<0||n<0)
		fault("mtx_resize:bad matrix dimension,m=%d,n=%d",m,n);
	if(r->lock&&(r->m!=m||r->n!=n))
		fault("mtx_resize:matrix dimensions locked",m,n);
#endif
	r->m=m;r->n=n;
	return r;
}

/*create new matrix object*/
mtx* mtx_create(int m,int n,int lock)
{
	mtx* p;
	p=malloc(sizeof(mtx));
	p->lock=0;
	mtx_resize(p,m,n);
	p->lock=lock;
	p->buf_len=0;
	p->data=NULL;
	p->t=MTX_NORMAL;
	p->pssm=NULL;
	return p;
}

/*print the matrix to stdout or FILE*/
void mtx_print(mtx* a)
{
	la_print(mtx_get_buf_read(a),a->m,a->n,NULL);
}
void mtx_fprint(FILE* fp,mtx* a)
{
	la_fprint(fp,mtx_get_buf_read(a),a->m,a->n,NULL);
}

/*r=a*/
mtx* mtx_copy(mtx* r,mtx* a)
{
	mtx_resize(r,a->m,a->n);
	if(a->m*a->n>0)
		la_copy(mtx_get_buf_write(r),
		mtx_get_buf_read(a),a->m,a->n);
	return r;
}

mtx* mtx_cover(mtx* r,mtx* a)
{
	int m,n,i,j;
#ifdef _SAFE_CHECK_
	if(mtx_is_valid(r)==0)
		fault("mtx_cover:ret matrix must be valid");
#endif
	m=r->m;if(m>a->m)m=a->m;
	n=r->n;if(n>a->n)n=a->n;
	for(i=1;i<=m;i++)
		for(j=1;j<=n;j++)
			mtx_set_element(r,i,j,mtx_get_element(a,i,j));
	return r;
}

mtx* mtx_copy_data(mtx* r,double* data,int m,int n)
{
	mtx_resize(r,m,n);
	la_copy(mtx_get_buf_write(r),
		data,m,n);
	return r;

}

/*r=diag(a,r)*/
mtx* mtx_diag(mtx* r,mtx* a,int k)
{
	int m,n;
	la_diag(NULL,a->data,a->m,a->n,k,&m,&n);
	mtx_resize(r,m,n);
	la_diag(mtx_get_buf_write(r),
		mtx_get_buf_read(a),a->m,a->n,k,NULL,NULL);
	return r;
}


/*r=a+b*/
mtx* mtx_add(mtx* r,mtx* a,mtx* b)
{
#ifdef _SAFE_CHECK_
	if(b->m!=a->m||b->n!=a->n)
		fault("mtx_add:Matrix dimensions must agree");
#endif
	mtx_resize(r,a->m,a->n);
	la_add(mtx_get_buf_write(r),
		mtx_get_buf_read(a),
		mtx_get_buf_read(b),a->m,a->n);
	return r;
}
/*a=a+b*/
mtx* mtx_add_self(mtx* a,mtx* b)
{
	return mtx_add(a,a,b);
}

/*r=a-b*/
mtx* mtx_minus(mtx* r,mtx* a,mtx* b)
{
#ifdef _SAFE_CHECK_
	if(b->m!=a->m||b->n!=a->n)
		fault("mtx_minus:Matrix dimensions must agree");
#endif
	mtx_resize(r,a->m,a->n);
	la_minus(mtx_get_buf_write(r),
		mtx_get_buf_read(a),
		mtx_get_buf_read(b),a->m,a->n);
	return r;
}
/*a=a-b*/
mtx* mtx_minus_self(mtx* a,mtx* b)
{
	return mtx_minus(a,a,b);
}
/*r=a.*b*/
mtx* mtx_dot_mul(mtx* r,mtx* a,mtx* b)
{
#ifdef _SAFE_CHECK_
	if(b->m!=a->m||b->n!=a->n)
		fault("mtx_dot_mul:Matrix dimensions must agree");
#endif
	mtx_resize(r,a->m,a->n);
	la_dot_multiply(mtx_get_buf_write(r),
		mtx_get_buf_read(a),
		mtx_get_buf_read(b),a->m,a->n);
	return r;
}
/*a=a.*b*/
mtx* mtx_dot_mul_self(mtx* a,mtx* b)
{
	return mtx_dot_mul(a,a,b);
}


/*r=a*b*/
mtx* mtx_mul(mtx* r,mtx* a,mtx* b)
{
#ifdef _SAFE_CHECK_
	if(a->n!=b->m)
		fault("矩阵不能相乘");
#endif
	mtx_resize(r,a->m,b->n);
	la_multiply(mtx_get_buf_write(r),
		mtx_get_buf_read(a),
		mtx_get_buf_read(b),a->m,a->n,b->n);
	return r;
}

mtx* mtx_mul_some(mtx* r,...)
{
	//mtx** ms;
	mtx* mt=(mtx*)1;
	va_list arg;
	fault("unsupport");
	va_start(arg,r);
	while(mt!=NULL)
	{
		mt=va_arg(arg,mtx*);
		printf("%X\n",mt);
	}
	va_end(arg);
	//free(ms);
	return r;
}

/*r=a*b*c*/
mtx* mtx_mul3(mtx* r,mtx* a,mtx* b,mtx* c)
{
#ifdef _SAFE_CHECK_
	if(a->n!=b->m||b->n!=c->m)
		fault("矩阵不能相乘");
#endif
	mtx_resize(r,a->m,c->n);
	la_multiply3(mtx_get_buf_write(r),
		mtx_get_buf_read(a),
		mtx_get_buf_read(b),
		mtx_get_buf_read(c),a->m,b->m,c->m,c->n);
	return r;
}

/*r=a*b*c*d*/
mtx* mtx_mul4(mtx* r,mtx* a,mtx* b,mtx* c,mtx* d)
{
#ifdef _SAFE_CHECK_
	if(a->n!=b->m||b->n!=c->m||c->n!=d->m)
		fault("矩阵不能相乘");
#endif
	mtx_resize(r,a->m,d->n);
	la_multiply4(mtx_get_buf_write(r),
		mtx_get_buf_read(a),
		mtx_get_buf_read(b),
		mtx_get_buf_read(c),
		mtx_get_buf_read(d),a->m,b->m,c->m,d->m,d->n);
	return r;
}
/*r=a*b*c*d*e*/
mtx* mtx_mul5(mtx* r,mtx* a,mtx* b,mtx* c,mtx* d,mtx* e)
{
#ifdef _SAFE_CHECK_
	if(a->n!=b->m||b->n!=c->m||c->n!=d->m||d->n!=e->m)
		fault("矩阵不能相乘");
#endif
	mtx_resize(r,a->m,e->n);
	la_multiply5(mtx_get_buf_write(r),
		mtx_get_buf_read(a),
		mtx_get_buf_read(b),
		mtx_get_buf_read(c),
		mtx_get_buf_read(d),
		mtx_get_buf_read(e),a->m,b->m,c->m,d->m,e->m,e->n);
	return r;
}
/*r=k*a*/
mtx* mtx_mul_number(mtx* r,double k,mtx* a)
{
	mtx_resize(r,a->m,a->n);
	la_multiply_number(mtx_get_buf_write(r),k,
		mtx_get_buf_read(a),a->m,a->n);
	return r;
}
/*a=k*a*/
mtx* mtx_mul_number_self(double k,mtx* a)
{
	return mtx_mul_number(a,k,a);
}

/*r=inv(a)*/
mtx* mtx_invert(mtx* r,mtx* a)
{
#ifdef _SAFE_CHECK_
	if(a->m!=a->n)
		fault("a不为方阵");
#endif
	mtx_resize(r,a->m,a->n);
	if(la_invert(mtx_get_buf_write(r),
		mtx_get_buf_read(a),a->m)==NULL)
		return NULL;/*矩阵不可逆*/
	return r;
}
mtx* mtx_invert_self(mtx* a)
{
#ifdef _SAFE_CHECK_
	if(a->m!=a->n)
		fault("a不为方阵");
#endif
	if(la_invert_self(mtx_get_buf_read(a),a->m)==NULL)
		return NULL;/*矩阵不可逆*/
	return a;
}
/*r=a'*/
mtx* mtx_transpose(mtx* r,mtx* a)
{
	mtx_resize(r,a->n,a->m);
	la_transpose(mtx_get_buf_write(r),
		mtx_get_buf_read(a),a->m,a->n);
	return r;
}

/*a=a'*/
mtx* mtx_transpose_self(mtx* a)
{
	int i;
	la_transpose_self(
		mtx_get_buf_read(a),a->m,a->n);
	i=a->m;
	a->m=a->n;
	a->n=i;
	return a;
}
double mtx_sum_element(mtx* p)
{
	int m=p->m,n=p->n,i,j;
	double* pd=mtx_get_buf_read(p);
	double sum=0;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			sum+=(*pd);
			pd++;
		}
	return sum;
}

mtx* mtx_sub(mtx* r,mtx* a,int m1,int m2,int n1,int n2)
{
	int m,n;
	la_sub(NULL,mtx_get_buf_read(a),a->m,a->n,m1,m2,n1,n2,&m,&n);
	if(m<1||n<1)return NULL;/*空*/
	mtx_resize(r,m,n);
	la_sub(mtx_get_buf_write(r),
		mtx_get_buf_read(a),a->m,a->n,
		m1,m2,n1,n2,NULL,NULL);
	return r;
}

mtx* mtx_tril(mtx* r,mtx* a,int k)
{
	mtx_resize(r,a->m,a->n);
	la_tril(mtx_get_buf_write(r),
		mtx_get_buf_read(a),a->m,a->n,k);
	return r;
}
mtx* mtx_tril_self(mtx* a,int k)
{
	la_tril_self(mtx_get_buf_read(a),a->m,a->n,k);
	return a;
}
mtx* mtx_triu(mtx* r,mtx* a,int k)
{
	mtx_resize(r,a->m,a->n);
	la_triu(mtx_get_buf_write(r),
		mtx_get_buf_read(a),a->m,a->n,k);
	return r;
}
mtx* mtx_triu_self(mtx* a,int k)
{
	la_triu_self(mtx_get_buf_read(a),a->m,a->n,k);
	return a;
}



/*初始化所有元素为0*/
mtx* mtx_ini_zeros(mtx* r,int m,int n)
{
	mtx_resize(r,m,n);
	la_zeros(mtx_get_buf_write(r),m,n);
	return r;
}

/*初始化矩阵为k倍的单位阵*/
mtx* mtx_ini_unit(mtx* r,int m,double k)
{
	mtx_resize(r,m,m);
	la_unit(mtx_get_buf_write(r),m,k);
	return r;
}

/*初始化所有元素为k*/
mtx* mtx_ini_ones(mtx* r,int m,int n,double k)
{
	mtx_resize(r,m,n);
	la_ones(mtx_get_buf_write(r),k,m,n);
	return r;
}


/*初始化所有元素为0-k的随机数*/
mtx* mtx_ini_rands(mtx* r,int m,int n,double k)
{
	mtx_resize(r,m,n);
	la_rands(mtx_get_buf_write(r),k,m,n);
	return r;
}

mtx* mtx_ini_diag(mtx* r,int m,double* data)
{
	int i;
	mtx_resize(r,m,m);
	la_zeros(mtx_get_buf_write(r),m,m);
	for(i=0;i<m;i++)
		mtx_set_element(r,i+1,i+1,data[i]);
	return r;
}
static void ORPlusConvolution( double *x, double *y, double *z, int Nx, int Ny, int *Nz )
/*  x: 0--Nx, y: 0--Ny,
    z: 0--Nz = Nx + Nz */
/*  z = x + y */
{
    register int i, j;
    int k;
	printf("<");
    *Nz = Nx + Ny;

    for ( i = 0; i <= *Nz; i ++ ) {
        z[i] = 0.0;
        for ( j = 0; j <= Ny; j ++ ) {
            k = i - j;
            if ( k < 0 ) break; /* because k is decreasing */
            else if ( k <= Nx &&x[k]!=0&&y[j]!=0) z[i] += x[k] * y[j];
//            if ( 0 <= k && k <= Nx ) z[i] += x[k] * y[j];
//        	fprintf( fpDebug, "i=%d  j=%d  k=%d\n  ", i, j, k );
            }
        }
	printf(">");

    return;
}

mtx* mtx_ORPlusConvolution(mtx* z,mtx* x,mtx* y)
{
	int nTemp;
	mtx_resize(z,x->m+y->m,1);
	ORPlusConvolution(mtx_get_buf_read(x),mtx_get_buf_read(y),mtx_get_buf_write(z),x->m,y->m,&nTemp);
	return z;
}

mtx* mtx_ORPlusConvolution_self(mtx* z,mtx* x)
{
	mtx* y=mtx_create(0,0,0);
	memcpy(y,z,sizeof(mtx));
	z->data=NULL;
	mtx_ORPlusConvolution(z,x,y);
	mtx_free(y);
	return z;
}

static int mtx_equal(double x,double y)
{
	return (fabs(x)-fabs(y))<1e-12;
}

double mtx_get_element(mtx* a,int m,int n)
{
#ifdef _SAFE_CHECK_
	if(m<1||n<1||m>a->m||n>a->n)
		fault("mtx_get_element:bad matrix dimension or exceed,m=%d,n=%d",m,n);
#endif
	return mtx_get_buf_read(a)[(m-1)*a->n+n-1];
}
void mtx_set_element(mtx* a,int m,int n,double e)
{
#ifdef _SAFE_CHECK_
	if(m<1||n<1||m>a->m||n>a->n)
		fault("mtx_get_element:bad matrix dimension or exceed,m=%d,n=%d",m,n);
#endif
	mtx_get_buf_write(a)[(m-1)*a->n+n-1]=e;
}

int mtx_get_nonzeros(mtx* a)
{
	int i,j,nz=0;
	if(a==NULL)return 0;
	if(a->t==MTX_SSM)return a->pssm->n;
	if(mtx_is_invalid(a))
		return 0;
	for(i=1;i<=a->m;i++)
		for(j=1;j<=a->n;j++)
			if(!mtx_equal(mtx_get_element(a,i,j),0))nz++;
	return nz;
}

int mtx_is_ssm(mtx* a)
{
	return a->t==MTX_SSM;
}

mtx* mtx_set_ssm(mtx* r,ssm* p)
{
	int i;
	i=0;
	mtx_clear_buf(r);
	r->t=MTX_SSM;
#ifdef _SAFE_CHECK_
	for(i=1;i<=p->n;i++)
		if (p->rn[i]>r->m || p->cn[i]>r->n)  fault("ERROR");
#endif
	r->pssm=p;
	return r;
}

ssm* mtx_convert_to_ssm(mtx* a)
{
	ssm* r;
	int i,j;
	double element;
	int nz=mtx_get_nonzeros(a);
	r=ssm_create(nz);
	for(i=1;i<=a->m;i++)
		for(j=1;j<=a->n;j++)
		{
			element=mtx_get_element(a,i,j);
			if(!mtx_equal(element,0))
				ssm_add(r,i,j,element);
		}
	mtx_clear_buf(a);
	a->t=MTX_SSM;
	a->pssm=r;
	return r;
}
/*------------------------------*------------------------------*/


/*create block matrix object*/
bmtx* bmtx_create(int m,int n,int rows_bks[],int cols_bks[])
{
	bmtx* p;
	int i;
#ifdef _SAFE_CHECK_
	if(m<1||n<1)
		fault("bmtx_new:bad matrix dimension,m=%d,n=%d",m,n);

	for(i=0;i<m;i++)
	{
		if(rows_bks[i]<0)
			fault("bmtx_new:bad block dimension,rows_bks[%d]=%d",i,rows_bks[i]);
	}
	for(i=0;i<n;i++)
	{
		if(cols_bks[i]<0)
			fault("bmtx_new:bad block dimension,cols_bks[%d]=%d",i,cols_bks[i]);
	}
#endif

	p=malloc(sizeof(bmtx));
	p->m=m;p->n=n;
	p->rm=0;p->rn=0;

	p->bks=malloc(sizeof(mtx*)*m*n);
	memset(p->bks,0,sizeof(mtx*)*m*n);
	p->bm=malloc(sizeof(int)*m);
	p->bn=malloc(sizeof(int)*n);
	p->sbm=malloc(sizeof(int)*m);
	p->sbn=malloc(sizeof(int)*n);

	memcpy(p->bm,rows_bks,sizeof(int)*m);
	memcpy(p->bn,cols_bks,sizeof(int)*n);

	for(i=0;i<m;i++)
	{
		p->sbm[i]=p->rm;
		p->rm+=rows_bks[i];
	}
	for(i=0;i<n;i++)
	{
		p->sbn[i]=p->rn;
		p->rn+=cols_bks[i];
	}

	return p;
}

/*get block*/
mtx* bmtx_get_bk_read(bmtx* p,int m,int n)
{
#ifdef _SAFE_CHECK_
	if(m<1||n<1||m>p->m||n>p->n)
		fault("bmtx_sub_mtx:bad matrix dimension or exceed,m=%d,n=%d",m,n);
#endif
	return p->bks[(m-1)*p->n+n-1];
}

mtx* bmtx_get_bk_write(bmtx* p,int m,int n)
{
	mtx** ppm;
#ifdef _SAFE_CHECK_
	if(m<1||n<1||m>p->m||n>p->n)
		fault("bmtx_sub_mtx:bad matrix dimension or exceed,m=%d,n=%d",m,n);
#endif
	ppm=&(p->bks[(m-1)*p->n+n-1]);
	if(*ppm==NULL)
		*ppm=mtx_create(bmtx_get_bk_rows(p,m),bmtx_get_bk_cols(p,n),1);
	return *ppm;
}

void bmtx_set_bk(bmtx* p,mtx* pm,int m,int n)
{
#ifdef _SAFE_CHECK_
	if(pm!=NULL
		&&bmtx_get_bk_rows(p,m)!=pm->m
		&&bmtx_get_bk_cols(p,n)!=pm->n)
		fault("ERROR");
	if(m<1||n<1||m>p->m||n>p->n)
		fault("bmtx_sub_mtx:bad matrix dimension or exceed,m=%d,n=%d",m,n);
#endif
	m--;n--;
	mtx_free(p->bks[m*p->n+n]);
	p->bks[m*p->n+n]=pm;
}

void bmtx_set_bk_ssm(bmtx* p,int m,int n,ssm* pssm)
{
	mtx_set_ssm(bmtx_get_bk_write(p,m,n),pssm);
}

void bmtx_clear_all(bmtx* p)
{
	int i,j;
	for(i=1;i<=p->m;i++)
		for(j=1;j<=p->n;j++)
			bmtx_clear_bk(p,i,j);
}

void bmtx_clear_bk(bmtx* p,int m,int n)
{
	bmtx_set_bk(p,NULL,m,n);
}


int bmtx_get_bk_rows(bmtx* p,int m)
{
	return p->bm[m-1];
}

int bmtx_get_bk_row_base(bmtx* p,int m)
{
	return p->sbm[m-1];
}


int bmtx_get_bk_cols(bmtx* p,int n)
{
	return p->bn[n-1];
}

int bmtx_get_bk_col_base(bmtx* p,int n)
{
	return p->sbn[n-1];
}

int bmtx_get_nonzeros(bmtx* p)
{
	int i,j,nz=0;
	for(i=1;i<=p->m;i++)
		for(j=1;j<=p->n;j++)
			nz+=mtx_get_nonzeros(bmtx_get_bk_read(p,i,j));
	return nz;
}

ssm* bmtx_create_ssm(bmtx* p)
{
	ssm* r;
	int i,j,k,l,idx=0,r_base,c_base;
	mtx* bk;ssm* bk_ssm;
	double d;
	int nz=bmtx_get_nonzeros(p);
	if(nz<1)return NULL;
	r=ssm_create(nz);
	for(i=1;i<=p->m;i++)
	{
		r_base=bmtx_get_bk_row_base(p,i);
		for(j=1;j<=p->n;j++)
		{
			c_base=bmtx_get_bk_col_base(p,j);
			bk=bmtx_get_bk_read(p,i,j);
			if(bk==NULL)continue;
			if(mtx_is_ssm(bk))
			{
				bk_ssm=bk->pssm;
				for(k=1;k<=bk_ssm->n;k++)
				{
					if(bk_ssm->a[k]!=0)
						ssm_add(r,r_base+bk_ssm->rn[k],
					c_base+bk_ssm->cn[k],bk_ssm->a[k]);
				}
				idx+=bk_ssm->n;
			}
			else if(mtx_is_valid(bk))
			{
				for(k=1;k<=bk->m;k++)
				{
					for(l=1;l<=bk->n;l++)
					{
						d=mtx_get_element(bk,k,l);
						if(!mtx_equal(d,0))
						{
							idx++;
#ifdef _SAFE_CHECK_
							if(idx>nz)
								fault("ERROR");
#endif
							if(d!=0)
								ssm_add(r,r_base+k,c_base+l,d);
						}
					}
				}
			}
		}
	}
#ifdef _SAFE_CHECK_
	if(idx!=nz)
		fault("ERROR");
#endif
	return r;
}

void bmtx_free(bmtx* p)
{
	int i,j;
	if(p==NULL)return;
	for(i=0;i<p->m;i++)
	{
		for(j=0;j<p->n;j++)
		{
			mtx_free(p->bks[i*p->n+j]);
		}
	}
	free(p->bks);
	free(p->bm);
	free(p->bn);
	free(p->sbm);
	free(p->sbn);
	free(p);
}


