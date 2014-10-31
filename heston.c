
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
double K=10;
void ConstructA1(int maxs,int maxv,double r,double **A1)
{	double k = 1.0/maxv;
	int i,j;
	for ( j=0;j<maxv*maxs;j++)
		for ( i = 0; i<maxs*maxv;i++)
		{
			A1[j][i]=0.0;
		}


	for (  j = 0;j<maxv;j++)
	{
		for( i =0;i<maxs;i++)
		{
			if(i==0 )
			{
				A1[j*maxs+i][j*maxs+i]= k*(j+1)*(i+1)*(i+1);
				A1[j*maxs+i][j*maxs+i+1] = -k*(j+1)*(i+1)*(i+1)/2-r*(i+1)/2;

			}
			else if (i==maxs-1 )
			{
				A1[j*maxs+i][j*maxs+i]= k*(j+1)*(i+1)*(i+1);
				A1[j*maxs+i][j*maxs+i-1] = -k*(j+1)*(i+1)*(i+1)/2;
			}
			else
			{
				A1[j*maxs+i][j*maxs+i]= k*(j+1)*(i+1)*(i+1);
				A1[j*maxs+i][j*maxs+i-1] = -k*(j+1)*(i+1)*(i+1)/2+r*(i+1)/2;
				A1[j*maxs+i][j*maxs+i+1] = -k*(j+1)*(i+1)*(i+1)/2-r*(i+1)/2;
			}

		}
	}

      //	return A1;
}
void ConstructA2(int maxs, int maxv,double **A2, double gama,double alfa,double beta,double k)
{

	int i,j;
	for (j=0;j<maxv*maxs;j++)
		for (  i = 0; i<maxs*maxv;i++)
		{
			A2[j][i]=0.0;
		}
	for ( j =0;j<maxv;j++)
		for( i = 0;i<maxs;i++)
		{

			A2[j*maxs+i][j*maxs+i]=gama*gama*(j+1)/k;
			if (j==0)
			{
				A2[(j)*maxs+i][(j+1)*maxs+i] = -gama*gama*(j+1)/2/k -alfa*beta/2/k+alfa*(j+1)/2;
			}
			else if(j==maxv-1)
			A2[(j)*maxs+i][(j-1)*maxs+i] =-gama*gama*(j+1)/k/2 ;
			else
			{
				A2[(j)*maxs+i][(j+1)*maxs+i] = -gama*gama*(j+1)/2/k -alfa*beta/2/k+alfa*(j+1)/2;
				A2[(j)*maxs+i][(j-1)*maxs+i] =-gama*gama*(j+1)/2/k +alfa*beta/2/k-alfa*(j+1)/2;
			}

		}
  printf("helloA2\n");

}
void ConstructA0(int maxs,int maxv,double **A0,double r ,double pho,double gama)
{

	int i,j;
	for (j=0;j<maxv*maxs;j++)
		for ( i = 0; i<maxs*maxv;i++)
		{
			A0[j][i]=0.0;
		}
	for ( j = 0;j<maxv;j++)
	{
		for ( i = 0;i<maxs;i++)
		{
			A0[j*maxs+i][j*maxs+i]=r;
			if(j==0 ){
				if (i==0)
				A0[j*maxs+i][(j+1)*maxs+i+1]=-pho*gama*(i+1)*(j+1)/4;
				else if (i==maxs-1)
				A0[j*maxs+i][(j+1)*maxs+i-1]=pho*gama*(i+1)*(j+1)/4;
				else
				{
					A0[j*maxs+i][(j+1)*maxs+i+1]=-pho*gama*(i+1)*(j+1)/4;
					A0[j*maxs+i][(j+1)*maxs+i-1]=pho*gama*(i+1)*(j+1)/4;
				}
			}
			else if(j==maxv-1 )
			{
				if (i==0)
				A0[j*maxs+i][(j-1)*maxs+i+1]=pho*gama*(i+1)*(j+1)/4;
				else if (i==maxs-1)
				A0[j*maxs+i][(j-1)*maxs+i-1]=-pho*gama*(i+1)*(j+1)/4;
				else
				{
					A0[j*maxs+i][(j-1)*maxs+i+1]=pho*gama*(i+1)*(j+1)/4;
					A0[j*maxs+i][(j-1)*maxs+i-1]=-pho*gama*(i+1)*(j+1)/4;
				}

			}
			else 
			{
				if (i==0){
				A0[j*maxs+i][(j-1)*maxs+i+1]=pho*gama*(i+1)*(j+1)/4;
				A0[j*maxs+i][(j+1)*maxs+i+1]=-pho*gama*(i+1)*(j+1)/4;
				}
				else if(i==maxs-1)
				{
					A0[j*maxs+i][(j-1)*maxs+i-1]=-pho*gama*(i+1)*(j+1)/4;
					A0[j*maxs+i][(j+1)*maxs+i-1]=pho*gama*(i+1)*(j+1)/4;
				}
				else
				{
				A0[j*maxs+i][(j+1)*maxs+i+1]=-pho*gama*(i+1)*(j+1)/4;
				A0[j*maxs+i][(j+1)*maxs+i-1]=pho*gama*(i+1)*(j+1)/4;
				A0[j*maxs+i][(j-1)*maxs+i-1]=-pho*gama*(i+1)*(j+1)/4;
				A0[j*maxs+i][(j-1)*maxs+i+1]=pho*gama*(i+1)*(j+1)/4;
				}
			}
		}
	}
	printf("helloA0\n");

}
double* Constructg(int maxs,double h)
{	int i;
	double * g=(double *)malloc(maxs*sizeof(double));
	for ( i = 0;i<maxs;i++){
		if ((K-(i+1)*h)>1e-7)
		{
			g[i]=K-(i+1)*h;
		}
		else g[i]=0.0;
		}
     //	printf("helloG!\n");
	return g;
}
double * ConstructU0(int maxs,int maxv,double *g)
{       int i,j;
	double * U0=(double *)malloc(maxs*maxv*sizeof(double)) ;
	for ( j=0;j<maxv;j++){
		for ( i = 0;i<maxs;i++)
		U0[(j*maxs+i)]=g[i];  }
	return U0;
}
void ConstructC11(int maxs,int maxv,double tao,double **A1,double ** A0,double **C11)
{
	int i,j;
	
	for ( i=0;i<maxv*maxs;i++)
		for ( j = 0; j<maxs*maxv;j++)
		{
			if (i==j)
			C11[i][j]=1.0+tao/2*(A1[i][j]+A0[i][j]/2.0);
			else
			C11[i][j]=tao/2*(A1[i][j]+A0[i][j]/2.0);
		}
	printf("helloC11\n");
	
}
double ** ConstructC12(int maxs,int maxv,double tao,double **A2,double ** A0)
{
	int i,j;
	double ** C12 =(double **) malloc(maxs*maxv*sizeof(double*));
	for( i =0;i<maxs*maxv;i++)
	C12[i]=(double *) malloc(maxs*maxv*sizeof(double));
	for (i=0;i<maxv*maxs;i++)
		for ( j = 0; j<maxs*maxv;j++)
		{
			if (i==j)
			C12[i][j]=1.0-tao/2*(A2[i][j]+A0[i][j]/2.0);
			else
			C12[i][j]=-tao/2*(A2[i][j]+A0[i][j]/2.0);

		}
	printf("helloC12\n");
	return C12;

}
double ** ConstructC21(int maxs,int maxv,double tao,double **A2,double ** A0)
{
	int i;
	int j;
	double ** C21 =(double **) malloc(maxs*maxv*sizeof(double*));
	for ( i = 0;i<maxs*maxv;i++)
	C21[i]=(double *) malloc(maxs*maxv*sizeof(double));
	for (i=0;i<maxv*maxs;i++)
		for ( j = 0; j<maxs*maxv;j++)
		{
			if (i==j)
			C21[i][j]=1.0+tao/2*(A2[i][j]+A0[i][j]/2.0);
			else
			C21[i][j]=tao/2*(A2[i][j]+A0[i][j]/2.0);
		}
	printf("helloC21\n");
	return C21;
}
double ** ConstructC22(int maxs,int maxv,double tao,double **A1,double ** A0)
{
	int i,j;
	double ** C22 = (double**)malloc(maxs*maxv*sizeof(double*));
	for( i =0;i<maxs*maxv;i++)
	C22[i]=(double *) malloc(maxs*maxv*sizeof(double));
	for ( i=0;i<maxv*maxs;i++)
		for ( j= 0; j<maxs*maxv;j++)
		{
			if (i==j)
			C22[i][j]=1.0-tao/2*(A1[i][j]+A0[i][j]/2.0);
			else
			C22[i][j]=-tao/2*(A1[i][j]+A0[i][j]/2.0);
		}
	printf("helloC22\n");
	return C22;
}
void  myUL(double ** C11,int maxs,int maxv,double **u,double**l)
{

	int m = maxs*maxv;
	int i,j,k;
    for ( i = 0;i<m;i++)//UL factorization
	{
	u[i][i]=1;
	l[m-1][i]=C11[m-1][i];
	u[i][m-1]=C11[i][m-1]/C11[m-1][m-1];
	}
    for ( i = m-2;i>=0;i--){
	for ( j =i;j>=0;j--){
	    double sum =0;
	    for ( k = i+1;k<m;k++)
		sum +=u[i][k]*l[k][j];

	    l[i][j]=C11[i][j]-sum;
	}
	for ( j = i;j>=0;j--){
	    double sum = 0;
	    for ( k = i+1;k<m;k++)
	    sum += u[j][k]*l[k][i];
	    u[j][i]= (C11[j][i]-sum)/l[i][i];
	}
	}
	printf("helloUL\n");

}
double ** inverseU(int maxs,int maxv,double **u)
{
	int m=maxs*maxv;
	int i,j,k;
       double **inverseu = (double**)malloc (m*sizeof(double*));
       double **bigU =(double **)malloc(m*sizeof(double*));
       for( i = 0;i<m;i++)
       {
		inverseu[i]=(double *) malloc(maxs*maxv*sizeof(double));
		bigU[i]=(double *) malloc(2*maxs*maxv*sizeof(double));
       }
	for ( i = 0;i<m;i++)
	{
	for ( j = 0;j<m;j++)
	{
		bigU[i][j]=u[i][j];
	}
	for (j = m;j<2*m;j++)
	{
		if(j-m==i)
		bigU[i][j]=1.0;
		else
		bigU[i][j]=0.0;
	}
	}
	for ( i =m-1;i>0;i--)
	{
	for ( j =i-1;j>=0;j--)//对角线上的元素
	{
		for (k = m+i;k<2*m;k++)//初等行变换
		{
		bigU[j][k]=bigU[j][k]-bigU[i][k]*bigU[j][i]/bigU[i][i];
		
		}
		bigU[j][i]=0.0;
	}
	}
	for( i =0;i<m;i++)
	{	for( j = 0;j<m;j++)
		inverseu[i][j]=bigU[i][j+m];
	}
	return inverseu;
}
void product(double **C,double *U,int maxs,int maxv, double *p)
{
	//double *D=(double *)malloc(maxs*maxv*sizeof(double));
	int m = maxs*maxv;
	int i,j;
	for (i = 0;i<m;i++){
		double sum = 0;
		for (j =0;j<m;j++)
		sum=sum+C[i][j]*U[j];
		p[i]=sum;
	}
}

double * solve(int maxs,int maxv,double *g,double **C12,double**C22,double **u1,double **u2,double ** l1,double **l2,int steps)
{
	int i,j,k,s;
	int m = maxs*maxv;
	double *U0=(double *)malloc(m*sizeof(double));
	double *U1=(double *)malloc(m*sizeof(double));
	
	U0=ConstructU0(maxs,maxv,g);
	double *q=(double *)malloc (m*sizeof(double));
	double *p = (double *)malloc (m*sizeof(double));
	
	double ** inverseu1;
	inverseu1=inverseU(maxs,maxv,u1);
	double ** inverseu2;
	inverseu2=inverseU(maxs,maxv,u2);
	for ( k =0;k<steps;k++)
	{
	product(C12,U0,maxs,maxv,q);
	product(inverseu1,q,maxs,maxv,p);
	for ( i =0;i<m;i++)
	{
		double sum = 0;
		for ( j = 0;j<i;j++)
		sum = sum+U1[j]*l1[i][j];
		U1[i]=(p[i]-sum)/l1[i][i];
		if (U1[i]-g[(i)%maxs]<1e-7)
		U1[i]=g[(i)%maxs];
		/*if(U1[i]<1e-6)
		{
			for (s=i;s<maxs*(i/maxs+1);s++)
			{
				U1[s]=0.0;
			}
			i=s;
			continue;
		}*/
	}
	
	for (i = 0;i<m;i++)
	U0[i]=U1[i];
	
	product(C22,U0,maxs,maxv,q);
	product(inverseu2,q,maxs,maxv,p);
	for ( i =0;i<m;i++)
	{
		double sum = 0;
		for ( j = 0;j<i;j++)
		sum = sum+U1[j]*l2[i][j];
		U1[i]=(p[i]-sum)/l2[i][i];
		if (U1[i]-g[(i)%maxs]<1e-7 )
		U1[i]=g[(i)%maxs];
		
		/*if(U1[i]<1e-6)
		{
			for (s=i;s<maxs*(i/maxs+1);s++)
			{
				U1[s]=0.0;
			}
			i=s;
			continue;
		}*/
	}
	for (i=0;i<m;i++)
	{
		
	}
	for (i = 0;i<m;i++)
	U0[i]=U1[i];
	
		printf("%d\n",k);
	}
	return U0;
}
int main()		
{
	int i,j,s;
	double solution;
	int maxs =80;int maxv = 64; int steps = 64;
	double t = 0.25;
	double h = 20.0/maxs;
	double k = 1.0/maxv;
	double S = 8.0;
    double	tao = 0.25/steps;
	double pho = 0.1;double gama = 0.9;double r = 0.1;
	double alfa = 5.0;double beta = 0.16;double sigma=0.0625;
	printf("hello1\n");
	double *g=Constructg(maxs,h);
	double **A1 =(double **) malloc(maxs*maxv*sizeof(double *));
	for ( i = 0;i<maxs*maxv;i++)
	{
		A1[i]=(double *)malloc(maxs*maxv*sizeof(double));
	}
	printf("A11\n");

	ConstructA1(maxs,maxv,r,A1);
	 printf("A1\n");

	double **A2 =(double **) malloc(maxs*maxv*sizeof(double *));
	for ( i = 0;i<maxs*maxv;i++)
	{
		A2[i]=(double *)malloc(maxs*maxv*sizeof(double));
	}
	ConstructA2(maxs,maxv,A2,gama,alfa,beta,k);
       double**	A0 =(double **) malloc(maxs*maxv*sizeof(double*));
	for ( i =0;i<maxs*maxv;i++)
	{
		A0[i]=(double *) malloc(maxs*maxv*sizeof(double));
	}
       ConstructA0(maxs,maxv,A0,r,pho,gama);
	   printf("====1 ===after construct a0\n");
	double **C11=(double **) malloc(maxs*maxv*sizeof(double *));
	printf("===== 2 ==after malloc c11 porinter a0\n");
	for ( i =0;i<maxs*maxv;i++)
	{
		C11[i]=(double *) malloc(maxs*maxv*sizeof(double));
	}

	ConstructC11(maxs,maxv,tao,A1,A0,C11);

	double **C12=(double **)ConstructC12(maxs,maxv,tao,A2,A0);

	double **C21=(double **)ConstructC21(maxs,maxv,tao,A2,A0);
	double **C22=(double **)ConstructC22(maxs,maxv,tao,A1,A0);
	
	free(A1);
	free(A2);
	free(A0);
	double **u1=(double **)malloc(maxs*maxv*sizeof(double *));
	double **l1=(double **)malloc(maxs*maxv*sizeof(double *));
	double **u2 =(double **) malloc(maxs*maxv*sizeof(double *));
	double **l2 =(double **) malloc(maxs*maxv*sizeof(double *));
	for ( i = 0;i<maxs*maxv;i++)
	{
		u1[i]=(double *) malloc(maxs*maxv*sizeof(double));
		u2[i]=(double *) malloc(maxs*maxv*sizeof(double));
		l1[i]=(double *) malloc(maxs*maxv*sizeof(double));
		l2[i]=(double *) malloc(maxs*maxv*sizeof(double));

	}
	
	myUL(C11,maxs,maxv,u1,l1);
	myUL(C21,maxs,maxv,u2,l2);
	free(C11);
	free(C21);
	double *myU0;//=(double *)malloc(maxs*maxv*sizeof(double));
	myU0=solve(maxs,maxv,g,C12,C22,u1,u2,l1,l2,steps);
	
	for( i =0;i<maxs;i++)
	{
		if (fabs(S-(i+1)*h)<h*0.25)
		break;
	}
	printf("i=%d\n",i);
	for (j = 0;j<maxv;j++)
	{
		if( fabs(sigma-(j+1)*k)<k*0.25 )
		break;
	}
	printf("j=%d\n",j);
	solution = myU0[j*maxs+i];
	for (s=0;s<maxs*maxv;s++){
	if(s%maxs==0)
	printf("%d\n%f\n ",s/maxs,myU0[s]);
	else
	printf("%f ",myU0[s]);
	}
	printf("solution is %f\n",solution);
	for (s=1;s<5;s++)
	printf("%f\n",myU0[j*maxs+i+4*s]);
	for (s=0;s<5;s++)
	printf("%f\n",myU0[(j+12)*maxs+i+4*s]);
	return 0;
      

}
