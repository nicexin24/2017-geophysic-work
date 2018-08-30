
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
#define PI 3.1415926
#define N 2000     //时间片数量
#define T 0.001    //每个时间片的时间间隔
#define t0 0.05    //时间延迟
#define X 350      //x单元数
#define Z 350      //z单元数
#define DX 10      //x轴步长
#define DZ 10      //z轴步长
#define M  5     //2M为空间差分阶数
#define L 10    //PML边界宽度
const int f0=30;  //主频
//定义d'函数
float dd(int x,int h,float v)
{
	if (h==1)
		return (L+M-x-1)*1.5*log(1000)*v/(pow(L,3))/DX/DX;
	else 
		return (x-X-L+M-1)*1.5*log(1000)*v/(pow(L,3))/DX/DX;
}
//定义d函数
float d(int x,int h,float v)
{
	if (h==1)
		return 1.5*log(1000)*v/(pow(L,3))*pow((L+M-x-1),2)/DX;
	else
		return 1.5*log(1000)*v/(pow(L,3))*pow((x-X-L+M-1),2)/DX;
}
//定义数组长度
const int X0=X+2*L;     
const int Z0=Z+2*L; 
int i,j,k;
//计算差分系数函数
void sc(double *c)   
{ 
	double a[M][M+1],l[M][M],t;
	int i,j,k,xj[M];
//计算增广矩阵
for (i = 0; i < M; i++) 
  { 
    for ( j = 0; j < M; j++) 
      a[i][ j] = pow( j + 1, 2 * (i + 1)); 
    if (i == 0) 
      a[i][M] = 1; 
    else 
      a[i][M] = 0; 
  } 
  for ( j = 0; j < M + 1; j++) 
    xj[j] = j; 
	for(k=0;k<M;k++)
	{
		if(a[k][k]==0)
			break;
		else
			for(i=k+1;i<M;i++)
			{
				l[i][k]=a[i][k]/a[k][k];
			for(j=k;j<=M;j++)
				a[i][j]=a[i][j]-l[i][k]*a[k][j];
			}
	}
	
	for(i=M-1;i>=0;i--)
	{	
		t=0;
		for(j=i+1;j<=M-1;j++)
		{	
			if(j==M)
			{
				t=0;
			}
			else
			t=t+a[i][j]*c[j];
		}
		c[i]=(a[i][M]-t)/a[i][i];
	}
	
}


//定义时间片
	  float p1[X0][Z0],p2[X0][Z0],p3[X0][Z0];   
	  float p11[X0][Z0],p12[X0][Z0],p13[X0][Z0]; 
	    float p22[X0][Z0],p23[X0][Z0];   
		 float p31[X0][Z0],p32[X0][Z0],p33[X0][Z0];   
		  float u1[X0][Z0],u2[X0][Z0],u3[X0][Z0];  
//定义pml计算函数
void pml(int i,int j,float a,float b,int x,int h,int dx ,int dz,float ss,float ff,float v){
					p13[i][j]=T*T*(v*v*a/(dx*dx)-d(x,h,v)*d(x,h,v)*p12[i][j]-2*d(x,h,v)*(p12[i][j]-p11[i][j])/T)+2*p12[i][j]-p11[i][j]+ss*ff;
					u3[i][j]=T*T*(-v*v*(p2[i][j]-p1[i][j])/dx*dd(x,h,v)-d(x,h,v)*d(x,h,v)*u2[i][j]-2*d(x,h,v)*(u2[i][j]-u1[i][j])/T)+2*u2[i][j]-u1[i][j]+ss*ff;
					p23[i][j]=(u3[i][j]*T+p22[i][j])/(1+T*d(x,h,v))+ss*ff;
					p33[i][j]=T*T*(v*v*b/(dz*dz))+2*p32[i][j]-p31[i][j]+ss*ff;
					p3[i][j]=p13[i][j]+p23[i][j]+p33[i][j]; 
}
int main()
{
cout<<"输入时间节点："<<endl;
int k0;
cin>>k0;
//计算具体差分系数
double c[M];
  sc(c); 
 //定义RICHER子波
	float f[N]; 
	int A=1;    //振幅
	
	for(k=0;k<N;k++)
	{
		f[k]=A*exp(-PI*PI*f0*f0*pow(k*T-t0,2))*(1-2*PI*PI*f0*f0*pow(k*T-t0,2));
	}
	
        //定义空间位置函数
	int s[X0][Z0]; 
	     for(i=0;i<X0;i++)
		for(j=0;j<Z0;j++)
		{
			if(i==X0/2 && j==Z0/2)
				s[i][j]=1;
			else
				s[i][j]=0;
		}

	
	//地下速度结构
	int V[X0][Z0];  
	for(i=0;i<X0;i++)
		for(j=0;j<Z0;j++)
		{
			V[i][j]=3000;
		}
		
		
      
//赋初值
	for(i=0;i<X0;i++)
		for(j=0;j<Z0;j++)
		{
			p2[i][j]=0;p3[i][j]=0;
			p12[i][j]=0;p13[i][j]=0;
                        p22[i][j]=0;p23[i][j]=0;
			p32[i][j]=0;p33[i][j]=0;
			u2[i][j]=0;u3[i][j]=0;
		}
	

//循环计算
	for(k=1;k<N;k++)         
	{
		 //时间片的递推 
        
				  for(i=0;i<X0;i++)          
				for(j=0;j<Z0;j++)
				{
					p1[i][j]=p2[i][j];
					p2[i][j]=p3[i][j];
					p11[i][j]=p12[i][j];
					p12[i][j]=p13[i][j];
					p22[i][j]=p23[i][j];
					p31[i][j]=p32[i][j];
					p32[i][j]=p33[i][j];	
					u1[i][j]=u2[i][j];
					u2[i][j]=u3[i][j];
				}
		//划分空间区域
		
				//正常区域
				
				for(i=M+L;i<L+X-M+1;i++)          
				      for(j=M+L;j<L+X-M+1;j++) 
				{
					float a(0),b(0);
					//计算空间导数算子
					for (int m=1;m<M+1;m++)
					{a=a+c[m-1]*(p2[i-m][j]-2*p2[i][j]+p2[i+m][j]);
						b=b+c[m-1]*(p2[i][j-m]-2*p2[i][j]+p2[i][j+m]);
					}	
	p3[i][j]=pow(V[i][j],2)*T*T*(a/(DX*DX)+b/(DZ*DZ))+2*p2[i][j]-p1[i][j]+s[i][j]*f[k];
				}
                
			//区域一 与x相关
			
					
				for(i=M;i<L+M;i++)          
				for(j=i;j<2*L+X-i;j++) 
				{
					float a(0),b(0);
					for (int m=1;m<M+1;m++)
					{a=a+c[m-1]*(p2[i-m][j]-2*p2[i][j]+p2[i+m][j]);
						b=b+c[m-1]*(p2[i][j-m]-2*p2[i][j]+p2[i][j+m]);
					}	
					
                                    pml(i,j,a,b,i,1,DX,DZ,s[i][j],f[k],V[i][j]);

				}
			
		   //区域二 与z相关
			
				  for(j=M;j<L+M;j++) 
					  for(i=j+1;i<Z+2*L-j+1;i++)          
				{
					float a(0),b(0);
					for (int m=1;m<M+1;m++)
					{
						a=a+c[m-1]*(p2[i-m][j]-2*p2[i][j]+p2[i+m][j]);
						b=b+c[m-1]*(p2[i][j-m]-2*p2[i][j]+p2[i][j+m]);
					}	
					 pml(i,j,b,a,j,1,DZ,DX,s[i][j],f[k],V[i][j]);


				}
			//区域三 与x相关
		
				for(i=L-M+X+1;i<2*L-M+X;i++)          
				for(j=X+2*L-i+1;j<i+1;j++) 
				{
					float a(0),b(0);
					for (int m=1;m<M+1;m++)
					{a=a+c[m-1]*(p2[i-m][j]-2*p2[i][j]+p2[i+m][j]);
						b=b+c[m-1]*(p2[i][j-m]-2*p2[i][j]+p2[i][j+m]);
					}	
					 pml(i,j,a,b,i,2,DX,DZ,s[i][j],f[k],V[i][j]);

				}
		  
			//区域四 与z相关
	
			
				  for(j=Z+L-M+1;j<Z+2*L-M;j++) 
					  for(i=2*L+X-j;i<j;i++)          
				{
					float a(0),b(0); 
					for (int m=1;m<M+1;m++)
					{a=a+c[m-1]*(p2[i-m][j]-2*p2[i][j]+p2[i+m][j]);
						b=b+c[m-1]*(p2[i][j-m]-2*p2[i][j]+p2[i][j+m]);
					}	
				pml(i,j,b,a,j,2,DZ,DX,s[i][j],f[k],V[i][j]);
				}
		
                  
		
	
              if(k==k0)
			
		{
			FILE *fp;
			if((fp=fopen("pml.dat","wb"))!=NULL)
			{
				for(i=0;i<X0;i++)
					for(j=0;j<Z0;j++)
						fwrite(&p3[i][j],sizeof(float),1,fp);
			}
                         else 
                      cout<<"fail to build new file";
			fclose(fp);exit(1);
		}
	}
	return 0;
}
