#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "la.h"


/*LINEAR ALGEBRA*/

/*
trnm

     Transpose a real square matrix in place A -> A~.

     void trnm(double *a,int n)
       a = pointer to array of n by n input matrix A
           This is overloaded by the transpose of A.
       n = dimension (dim(a)=n*n)
*/
static void trnm(double *a,int n)
{ double s,*p,*q;
  int i,j,e;
  for(i=0,e=n-1; i<n-1 ;++i,--e,a+=n+1){
    for(p=a+1,q=a+n,j=0; j<e ;++j){
      s= *p; *p++ = *q; *q=s; q+=n;
     }
   }
}

/*
mattr

     Transpose an m by n matrix A = B~.

     void mattr(double *a,double *b,int m,int n)
       a = pointer to array containing output n by m matrix 
       b = pointer to array containing input m by n matrix
            (matrices stored in row order)
       m,n = dimension parameters (dim(a)=dim(b)=n*m)
*/
static void mattr(double *a,double *b,int m,int n)
{ double *p; int i,j;
  for(i=0; i<n ;++i,++b)
    for(j=0,p=b; j<m ;++j,p+=n) *a++ = *p;
}



/*
mmul

     Multiply two real square matrices C = A * B.

     void mmul(double *c,double *a,double *b,int n)
     double *a,*b,*c; int n;
       a = pointer to store for left product matrix
       b = pointer to store for right product matrix
       c = pointer to store for output matrix
       n = dimension (dim(a)=dim(b)=dim(c)=n*n)
*/
static void mmul(double *c,double *a,double *b,int n)
{ double *p,*q,s; int i,j,k;
  trnm(b,n);
  for(i=0; i<n ;++i,a+=n){
    for(j=0,q=b; j<n ;++j){
      for(k=0,p=a,s=0.; k<n ;++k) s+= *p++ * *q++;
      *c++ =s;
     }
   }
  trnm(b,n);
}



/*
rmmult

     Multiply two matrices Mat = A*B.

     void rmmult(double *mat,double *a,double *b,int m,int k,int n)
     double mat[],a[],b[]; int m,k,n;
       mat = array containing m by n product matrix at exit
       a = input array containing m by k matrix
       b = input array containing k by n matrix
            (all matrices stored in row order)
       m,k,n = dimension parameters of arrays
*/
static void rmmult(double *rm,double *a,double *b,int n,int m,int l)
{ double z,*q0,*p,*q; int i,j,k;
  q0=(double *)calloc(m,sizeof(double));
  for(i=0; i<l ;++i,++rm){
    for(k=0,p=b+i; k<m ;p+=l) q0[k++]= *p;
    for(j=0,p=a,q=rm; j<n ;++j,q+=l){
      for(k=0,z=0.; k<m ;) z+= *p++ * q0[k++];
      *q=z;
     }
   }
  free(q0);
}


/*
vmul

     Multiply a vector by a matrix Vp = Mat*V.

     void vmul(double *vp,double *mat,double *v,int n)
       vp = pointer to array containing output vector
       mat = pointer to array containing input matrix in row order
       v = pointer to array containing input vector
       n = dimension of vectors (mat is n by n)
*/
static void vmul(double *vp,double *mat,double *v,int n)
{ double s,*q; int k,i;
  for(k=0; k<n ;++k){
    for(i=0,q=v,s=0.; i<n ;++i) s+= *mat++ * *q++;
    *vp++ =s;
   }
}


/*
vnrm

     Compute the inner product of two real vectors, p = u~*v.

     double vnrm(double *u,double *v,int n)
       u = pointer to array of input vector u
       v = pointer to array of input vector v
       n = dimension (dim(u)=dim(v)=n)
      return value: p = u~*v (dot product of u and v)
*/
static double vnrm(double *u,double *v,int n)
{ double s; int i;
  for(i=0,s=0.; i<n ;++i) s+= *u++ * *v++;
  return s;
}


/*
otrma

     Perform an orthogonal similarity transform C = A*B*A~.

     void otrma(double *c,double *a,double *b,int n)
       c = pointer to array of output matrix C
       a = pointer to array of transformation A
       b = pointer to array of input matrix B
       n = dimension (dim(a)=dim(b)=dim(c)=n*n)
*/
static void otrma(double *c,double *a,double *b,int n)
{ double z,*q0,*p,*s,*t;
  int i,j,k;
  q0=(double *)calloc(n,sizeof(double));
  for(i=0; i<n ;++i,++c){
    for(j=0,t=b; j<n ;++j){
      for(k=0,s=a+i*n,z=0.; k<n ;++k) z+= *t++ * *s++;
      q0[j]=z;
     }
    for(j=0,p=c,t=a; j<n ;++j,p+=n){
      for(k=0,s=q0,z=0.; k<n ;++k) z+= *t++ * *s++;
      *p=z;
     }
   }
  free(q0);
}

/*
otrsm

     Perform a similarity transform on a symmetric matrix S = A*B*A~.

     void otrsm(double *sm,double *a,double *b,int n)
       sm = pointer to array of output matrix S
       a = pointer to array of transformation matrix A
       b = pointer to array of symmetric input matrix B
       n = dimension (dim(a)=dim(b)=dim(sm)=n*n)
*/
static void otrsm(double *sm,double *a,double *b,int n)
{ double z,*q0,*p,*s,*t;
  int i,j,k;
  q0=(double *)calloc(n,sizeof(double));
  for(i=0; i<n ;++i){
    for(j=0,t=b; j<n ;++j){
      for(k=0,s=a+i*n,z=0.; k<n ;++k) z+= *t++ * *s++;
      q0[j]=z;
     }
    for(j=0,p=sm+i,t=a; j<=i ;++j,p+=n){
      for(k=0,s=q0,z=0.; k<n ;++k) z+= *t++ * *s++;
      *p=z; if(j<i) sm[i*n+j]=z;
     }
   }
  free(q0);
}


/*
minv

     Invert (in place) a general real matrix A -> Inv(A).

     int minv(double a[],int n)
       a = array containing the input matrix A
           This is converted to the inverse matrix.
       n = dimension of the system (i.e. A is n x n )
      return: 0 -> normal exit
              1 -> singular input matrix
*/
static int minv(double *a,int n)
{ int lc,*le; double s,t,tq=0.,zr=1.e-15;
  double *pa,*pd,*ps,*p,*q,*q0;
  int i,j,k,m;
  le=(int *)malloc(n*sizeof(int));
  q0=(double *)malloc(n*sizeof(double));
  for(j=0,pa=pd=a; j<n ;++j,++pa,pd+=n+1){
    if(j>0){
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
      for(i=1; i<n ;++i){ lc=i<j?i:j;
        for(k=0,p=pa+i*n-j,q=q0,t=0.; k<lc ;++k) t+= *p++ * *q++;
      	q0[i]-=t;
       }
      for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
     }
    s=fabs(*pd); lc=j;
    for(k=j+1,ps=pd; k<n ;++k){
      if((t=fabs(*(ps+=n)))>s){ s=t; lc=k;}
     }
    tq=tq>s?tq:s; if(s<zr*tq){ free(le-j); free(q0); return -1;}
    *le++ =lc;
    if(lc!=j){
      for(k=0,p=a+n*j,q=a+n*lc; k<n ;++k){
        t= *p; *p++ = *q; *q++ =t;
       }
     }
    for(k=j+1,ps=pd,t=1./ *pd; k<n ;++k) *(ps+=n)*=t;
    *pd=t;
   }
  for(j=1,pd=ps=a; j<n ;++j){
    for(k=0,pd+=n+1,q= ++ps; k<j ;++k,q+=n) *q*= *pd;
   }
  for(j=1,pa=a; j<n ;++j){ ++pa;
    for(i=0,q=q0,p=pa; i<j ;++i,p+=n) *q++ = *p;
    for(k=0; k<j ;++k){ t=0.;
      for(i=k,p=pa+k*n+k-j,q=q0+k; i<j ;++i) t-= *p++ * *q++;
      q0[k]=t;
     }
    for(i=0,q=q0,p=pa; i<j ;++i,p+=n) *p= *q++;
   }
  for(j=n-2,pd=pa=a+n*n-1; j>=0 ;--j){ --pa; pd-=n+1;
    for(i=0,m=n-j-1,q=q0,p=pd+n; i<m ;++i,p+=n) *q++ = *p;
    for(k=n-1,ps=pa; k>j ;--k,ps-=n){ t= -(*ps);
      for(i=j+1,p=ps,q=q0; i<k ;++i) t-= *++p * *q++;
      q0[--m]=t;
     }
    for(i=0,m=n-j-1,q=q0,p=pd+n; i<m ;++i,p+=n) *p= *q++;
   }
  for(k=0,pa=a; k<n-1 ;++k,++pa){
    for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *q++ = *p;
    for(j=0,ps=a; j<n ;++j,ps+=n){
      if(j>k){ t=0.; p=ps+j; i=j;}
      else{ t=q0[j]; p=ps+k+1; i=k+1;}
      for(; i<n ;) t+= *p++ *q0[i++];
      q0[j]=t;
     }
    for(i=0,q=q0,p=pa; i<n ;++i,p+=n) *p= *q++;
   }
  for(j=n-2,le--; j>=0 ;--j){
    for(k=0,p=a+j,q=a+ *(--le); k<n ;++k,p+=n,q+=n){
      t=*p; *p=*q; *q=t;
     }
   }
  free(le); free(q0);
  return 0;
}


/*
psinv

     Invert (in place) a symmetric real matrix, V -> Inv(V).

     int psinv(double v[],int n)
       v = array containing a symmetric input matrix
           This is converted to the inverse matrix.
       n = dimension of the system (dim(v)=n*n)
      return: 0 -> normal exit
              1 -> input matrix not positive definite

           The input matrix V is symmetric (V[i,j] = V[j,i]).
*/
static int psinv(double *v,int n)
{ double z,*p,*q,*r,*s,*t; int j,k;
  for(j=0,p=v; j<n ;++j,p+=n+1){
    for(q=v+j*n; q<p ;++q) *p-= *q* *q;
    if(*p<=0.) return -1;
    *p=sqrt(*p);
    for(k=j+1,q=p+n; k<n ;++k,q+=n){
      for(r=v+j*n,s=v+k*n,z=0.; r<p ;) z+= *r++ * *s++;
      *q-=z; *q/= *p;
     }
   }
  trnm(v,n);
  for(j=0,p=v; j<n ;++j,p+=n+1){ *p=1./ *p;
    for(q=v+j,t=v; q<p ;t+=n+1,q+=n){
      for(s=q,r=t,z=0.; s<p ;s+=n) z-= *s * *r++;
      *q=z* *p; }
   }
  for(j=0,p=v; j<n ;++j,p+=n+1){
    for(q=v+j,t=p-j; q<=p ;q+=n){
      for(k=j,r=p,s=q,z=0.; k<n ;++k) z += *r++ * *s++;
      *t++ =(*q=z); }
   }
  return 0;
}




/*
ruinv

     Invert an upper right triangular matrix T -> Inv(T).

     int ruinv(double *a,int n)
       a = pointer to array of upper right triangular matrix
           This is replaced by the inverse matrix.
       n = dimension (dim(a)=n*n)
      return value: status flag, with 0 -> matrix inverted
                                     -1 -> matrix singular
*/
static int ruinv(double *a,int n)
{ int j; double fabs();
  double tt,z,*p,*q,*r,*s,*t;
  for(j=0,tt=0.,p=a; j<n ;++j,p+=n+1) if((z=fabs(*p))>tt) tt=z;
  tt*=1.e-16;
  for(j=0,p=a; j<n ;++j,p+=n+1){
    if(fabs(*p)<tt) return -1;
    *p=1./ *p;
    for(q=a+j,t=a; q<p ;t+=n+1,q+=n){
      for(s=q,r=t,z=0.; s<p ;s+=n) z-= *s * *r++;
      *q=z* *p;
     }
   }
  return 0;
}






double* la_zeros(double* a,int m,int n)
{
	memset(a,0,sizeof(double)*m*n);
	return a;
}

double* la_ones(double *a,double e,int m,int n)
{
	int i;
	double *p=a;
	for(i=0;i<m*n;i++)
	{
		*p=e;
		p++;
	}
	return a;
}

double* la_unit(double *a,int n,double e)
{
	int i;
	la_zeros(a,n,n);
	for(i=0;i<n;i++)
	{
		*(a+i*n+i)=e;
	}
	return a;
}

double* la_rands(double *a,double k,int m,int n)
{
	int i;
	double *p=a;
	srand( (unsigned)time( NULL ) );
	for(i=0;i<m*n;i++)
	{
		*p=k*rand()/RAND_MAX;
		p++;
	}
	return a;
	
}

double* la_sub(double* b,double* a,int m,int n,int m1,int m2,int n1,int n2,int* mb,int* nb)
{
	int s,ns,i;
	if(m1<1)m1=1;
	if(m2>m)m2=m;
	if(n1<1)n1=1;
	if(n2>n)n2=n;
	m1--;m2--;n1--;n2--;
	ns=n2-n1+1;
	if(b!=NULL)
	{
		s=sizeof(double)*ns;
		for(i=m1;i<=m2;i++)
		{
			memcpy(b+(i-m1)*ns,a+i*n+n1,s);
		}
	}
	if(mb!=NULL)*mb=(m2-m1+1);
	if(nb!=NULL)*nb=ns;
	return b;
}
/*b可以为NULL，此时可以获取得到的数组大小*/
double* la_diag(double* b ,double* a,int m,int n,int r,int* pmb,int* pnb)
{
	/*
	n-r>0
	m+r>0
	*/
	int mb,nb,i;
	if(m==1||n==1)
	{
		mb=m+n-1+abs(r);
		nb=mb;
		if(b!=NULL)
		{
			la_zeros(b,mb,nb);
			for(i=0;i<m+n-1;i++)
			{
				if(r>=0)*(b+i*nb+i+r)=*(a+i);
				else *(b+(i-r)*nb+i)=*(a+i);
			}
		}
	}
	else
	{
		nb=1;
		if(r>=0)
		{
			mb=n-r;
			if(mb>m)mb=m;
		}
		else
		{
			mb=m+r;
			if(mb>n)mb=n;
		}
		if(b!=NULL)
		{
			for(i=0;i<mb;i++)
			{
				if(r>=0)*(b+i)=*(a+i*n+i+r);
				else *(b+i)=*(a+(i-r)*n+i);
			}
		}
	}
	if(pmb!=NULL)*pmb=mb;
	if(pnb!=NULL)*pnb=nb;
	return b;
}


double* la_triu_self(double *a,int m,int n,int r)
{
	int i,k;
	if(r>=0)
	{
		for(i=0;i<m;i++)
		{
			k=r+i;
			if(k>n)k=n;
			memset(a+i*n,0,sizeof(double)*k);
		}
	}
	else
	{
		r=-r;
		for(i=r+1;i<m;i++)
		{
			memset(a+i*n,0,sizeof(double)*(m-r-1));
		}
	}
	return a;
}
double* la_triu(double* b,double *a,int m,int n,int r)
{
	la_copy(b,a,m,n);
	la_triu_self(b,m,n,r);
	return b;
}

double* la_tril_self(double *a,int m,int n,int r)
{
	int i,k;
	if(r>=0)
	{
		k=n;
		if(k>m+r+1)k=m+r+1;
		for(i=r+1;i<k;i++)
		{
			memset(a+(i-r-1)*n+i,0,sizeof(double)*(n-i));
		}
	}
	else
	{
		r=-r;
		k=m;
		if(k>n+r)k=n+r;
		for(i=0;i<k;i++)
		{
			memset(a+i*n+(i>=r?(i-r+1):0),0,sizeof(double)*(i>=r?(n+r-i-1):n));
		}
	}
	return a;
}
double* la_tril(double* b,double *a,int m,int n,int r)
{
	la_copy(b,a,m,n);
	la_tril_self(b,m,n,r);
	return b;
}


void la_print(double *a,int m,int n,char *fmt)
{
	la_fprint(stdout,a,m,n,fmt);
}
void la_fprint(void *fp,double *a,int m,int n,char *fmt)
{
	int i,j; 
	double *p;
	if(fmt==NULL)fmt="%8.4g";
	for(i=0,p=a; i<m ;++i)
	{
		for(j=0; j<n ;++j)
		{
			fprintf((FILE*)fp,fmt,*p++);
			if(j<n-1)fprintf((FILE*)fp,",");
		}
		/*if(n>1)fprintf((FILE*)fp,"\n");
		elsefprintf((FILE*)fp,";");*/
		fprintf((FILE*)fp,"\n");
	}
	//fprintf((FILE*)fp,"\n");
}





/*b=a*/
double* la_copy(double* b,const double* a,int m,int n)
{
	memcpy(b,a,sizeof(double)*m*n);
	return b;
}
/*c=a*b*/
double* la_multiply(double *c,double *a,double *b,int m,int k,int n)
{
	if(m==1&&n==1)
	{
		*c=vnrm(a,b,k);
	}
	else if(m==k&&n==1)
	{
		vmul(c,a,b,m);
	}
	else if(m==k&&k==n)
	{
		mmul(c,a,b,n);
	}
	else rmmult(c,a,b,m,k,n);
	return c;
}
double* la_multiply3(double *d,double *a,double *b,double *c,int m,int k,int l,int n)
{
	double* tmp;
	if(m*l<=k*n)
	{
		tmp=malloc(sizeof(double)*m*l);
		la_multiply(tmp,a,b,m,k,l);
		la_multiply(d,tmp,c,m,l,n);
	}
	else
	{
		tmp=malloc(sizeof(double)*k*n);
		la_multiply(tmp,b,c,k,l,n);
		la_multiply(d,a,tmp,m,k,n);
	}
	free(tmp);
	return d;
}

double* la_multiply4(double *e,double *a,double *b,double *c,double *d,int m,int k,int l,int x,int n)
{
	double* tmp=malloc(sizeof(double)*m*x);
	la_multiply3(tmp,a,b,c,m,k,l,x);
	la_multiply(e,tmp,d,m,x,n);
	free(tmp);
	return e;
}

double* la_multiply5(double *f,double *a,double *b,double *c,double *d,double *e
					 ,int m,int k,int l,int x,int y,int n)
{
	double* tmp=malloc(sizeof(double)*m*x);
	la_multiply3(tmp,a,b,c,m,k,l,x);
	la_multiply3(f,tmp,d,e,m,x,y,n);
	free(tmp);
	return f;
}



/*c=a.b*/
double* la_dot_multiply(double *c,double *a,double *b,int m,int n)
{
	int i;
	double *pa=a,*pb=b,*pc=c;
	for(i=0;i<m*n;i++)
	{
		*pc=(*pa)*(*pb);
		pa++;pb++;pc++;
	}
	return c;
}
/*a=a.b*/
double* la_dot_multiply_self(double *a,double *b,int m,int n)
{
	return la_dot_multiply(a,a,b,m,n);
}



/*c=a+b*/
double* la_add(double* c,double* a,double* b,int m,int n)
{
	int i;
	double *pa=a,*pb=b,*pc=c;
	for(i=0;i<m*n;i++)
	{
		*pc=(*pa)+(*pb);
		pa++;pb++;pc++;
	}
	return c;
}
/*a=a+b*/
double* la_add_self(double* a,double*b,int m,int n)
{
	return la_add(a,a,b,m,n);
}

/*c=a-b*/
double* la_minus(double* c,double* a,double* b,int m,int n)
{
	int i;
	double *pa=a,*pb=b,*pc=c;
	for(i=0;i<m*n;i++)
	{
		*pc=(*pa)-(*pb);
		pa++;pb++;pc++;
	}
	return c;
}
/*a=a-b*/
double* la_minus_self(double* a,double*b,int m,int n)
{
	return la_minus(a,a,b,m,n);
}


/*b=k*a*/
double* la_multiply_number(double* c,double k,double* a,int m,int n)
{
	int i;
	double *pa=a,*pc=c;
	for(i=0;i<m*n;i++)
	{
		*pc=k*(*pa);
		pa++;pc++;
	}
	return c;
}

/*a=k*a*/
double* la_multiply_number_self(double k,double* a,int m,int n)
{
	return la_multiply_number(a,k,a,m,n);
}

/*a=inv(a)*/
double* la_invert_self(double* a,int n)
{
	if(minv(a,n)!=0)return NULL;
	else return a;
}
/*b=inv(a)*/
double* la_invert(double* b,double* a,int n)
{
	la_copy(b,a,n,n);
	return la_invert_self(b,n);
}

/*a=a'*/
double* la_transpose_self(double *a,int m,int n)
{
	double* tmp;
	if(m==n)
	{
		trnm(a,n);
		return a;
	}
	tmp=malloc(sizeof(double)*m*n);
	mattr(tmp,a,m,n);
	la_copy(a,tmp,m,n);
	free(tmp);
	return a;
}
/*b=a'*/
double* la_transpose(double* b,double *a,int m,int n)
{
	if(m==n)
	{
		la_copy(b,a,n,n);
		trnm(b,n);
		return b;
	}
	mattr(b,a,m,n);
	return b;
}

#if 0

/*test*/
void la_main(void)
{
	int m=3,n=4,k=5,mb,nb;
	double *a,*b,*c,e=2.7;

	a=malloc(sizeof(double)*m*n);
	b=malloc(sizeof(double)*n*k);
	c=malloc(sizeof(double)*m*k);




	printf("分配%d×%d的各个元素为零矩阵\n",m,n);
	la_print(la_zeros(a,m,n),m,n,NULL);
	printf("\n\n分配%d×%d的各个元素为%f矩阵\n",m,n,e);
	la_print(la_ones(a,e,m,n),m,n,NULL);
	printf("\n\n分配%d×%d的对角元素为%f矩阵\n",m,m,e);
	la_print(la_unit(a,m,e),m,m,NULL);
	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",m,n);
	la_print(la_rands(a,e,m,n),m,n,NULL);

	printf("\n\n获取第2-3列\n",m,n);
	la_print(la_sub(b,a,m,n,0,0,2,3,NULL,NULL),m,2,NULL);
	printf("\n\n获取第2-3行\n",m,n);
	la_print(la_sub(b,a,m,n,2,3,0,0,NULL,NULL),2,n,NULL);

	printf("\n\n获取第-1对角\n",m,n);
	la_diag(b,a,m,n,-1,&mb,&nb);
	la_print(b,mb,nb,NULL);

	printf("\n\n获取第0对角\n",m,n);
	la_diag(b,a,m,n,0,&mb,&nb);
	la_print(b,mb,nb,NULL);

	printf("\n\n获取第1对角\n",m,n);
	la_diag(b,a,m,n,1,&mb,&nb);
	la_print(b,mb,nb,NULL);

	printf("\n\n获取第2对角\n",m,n);
	la_diag(b,a,m,n,2,&mb,&nb);
	la_print(b,mb,nb,NULL);

	printf("\n\n向量\n",m,n);
	la_print(la_rands(a,e,m,1),m,1,NULL);
	
	printf("\n\n生产-1对角阵\n",m,n);
	la_diag(b,a,m,1,-1,&mb,&nb);
	la_print(b,mb,nb,NULL);

	printf("\n\n生产0对角阵\n",m,n);
	la_diag(b,a,m,1,0,&mb,&nb);
	la_print(b,mb,nb,NULL);

	printf("\n\n生产1对角阵\n",m,n);
	la_diag(b,a,m,1,1,&mb,&nb);
	la_print(b,mb,nb,NULL);

	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",m,n);
	la_print(la_rands(a,e,m,n),m,n,NULL);
	printf("\n\nC=3*A\n",m,n);
	la_print(la_multiply_number(c,3,a,m,n),m,n,NULL);
	printf("\n\nA=3*A\n",m,n);
	la_print(la_multiply_number_self(3,a,m,n),m,n,NULL);

	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",m,n);
	la_print(la_rands(b,e,m,n),m,n,NULL);
	printf("\n\nC=A+B\n",m,n);
	la_print(la_add(c,a,b,m,n),m,n,NULL);
	printf("\n\nA=A+B\n",m,n);
	la_print(la_add_self(a,b,m,n),m,n,NULL);
	printf("\n\nC=A-B\n",m,n);
	la_print(la_minus(c,a,b,m,n),m,n,NULL);
	printf("\n\nA=A-B\n",m,n);
	la_print(la_minus_self(a,b,m,n),m,n,NULL);
	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",m,n);
	la_print(la_rands(a,e,m,n),m,n,NULL);
	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",n,k);
	la_print(la_rands(b,e,n,k),n,k,NULL);
	printf("\n\nC=A*B\n",m,n);
	la_print(la_multiply(c,a,b,m,n,k),m,k,NULL);

	
	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",m,n);
	la_print(la_rands(a,e,m,n),m,n,NULL);

	printf("\n\nA=A'\n",m,n);
	la_print(la_transpose_self(a,m,n),n,m,NULL);	
	
	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",m,n);
	la_print(la_rands(a,e,m,n),m,n,NULL);
	printf("\n\nA=la_triu(3)\n",m,n);
	la_print(la_triu(b,a,m,n,3),m,n,NULL);	
	printf("\n\nA=la_triu(2)\n",m,n);
	la_print(la_triu(b,a,m,n,2),m,n,NULL);	
	printf("\n\nA=la_triu(1)\n",m,n);
	la_print(la_triu(b,a,m,n,1),m,n,NULL);	
	printf("\n\nA=la_triu(0)\n",m,n);
	la_print(la_triu(b,a,m,n,0),m,n,NULL);	
	printf("\n\nA=la_triu(-1)\n",m,n);
	la_print(la_triu(b,a,m,n,-1),m,n,NULL);	
	printf("\n\nA=la_triu(-2)\n",m,n);
	la_print(la_triu(b,a,m,n,-2),m,n,NULL);	

	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",m,n);
	la_print(la_rands(a,e,m,n),m,n,NULL);
	printf("\n\nA=la_tril(3)\n",m,n);
	la_print(la_tril(b,a,m,n,3),m,n,NULL);	
	printf("\n\nA=la_tril(2)\n",m,n);
	la_print(la_tril(b,a,m,n,2),m,n,NULL);	
	printf("\n\nA=la_tril(1)\n",m,n);
	la_print(la_tril(b,a,m,n,1),m,n,NULL);	
	printf("\n\nA=la_tril(0)\n",m,n);
	la_print(la_tril(b,a,m,n,0),m,n,NULL);	
	printf("\n\nA=la_tril(-1)\n",m,n);
	la_print(la_tril(b,a,m,n,-1),m,n,NULL);	
	printf("\n\nA=la_tril(-2)\n",m,n);
	la_print(la_tril(b,a,m,n,-2),m,n,NULL);	


	printf("\n\n分配%d×%d的各个元素为随机数矩阵\n",m,m);
	la_print(la_rands(a,e,m,n),m,m,NULL);
	printf("\n\nC=INV(A)\n",m,m);
	la_invert(c,a,m);
	la_print(c,m,m,NULL);
	printf("\n\nC=A*INV(A)\n",m,m);
	la_print(la_multiply(b,c,a,m,m,m),m,m,NULL);
	printf("\n\nA=INV(A)\n",m,m);
	la_invert_self(a,m);
	la_print(a,m,m,NULL);
	printf("\n\nA=A'\n",m,m);
	la_print(la_transpose_self(a,m,m),m,m,NULL);	

	free(a);
	free(b);
	free(c);
}

#endif