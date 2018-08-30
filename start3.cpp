#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
#define PI 3.141593
#define N 2000     //时间片数量
#define T 0.001    //每个时间片的时间间隔
#define X 500      //x单元数
#define Z 500      //z单元数
#define DX 10      //x轴步长
#define DZ 10      //z轴步长
#define XH 10      //x轴外延
#define ZH 10      //z轴外延


int main()
{
	int X0=X+XH;     //数组总长度
	int Z0=Z+ZH;     //数组总长度
	float f[N]; //RICHER子波
	int f0=30;  //主频
	int A=1;    //振幅
	int t0=0.03;//延迟时间
	int i,j,k;
	for(k=0;k<N;k++)
	{
		f[k]=A*exp(-PI*PI*f0*f0*(k*T-t0)*(k*T-t0))*(1-2*PI*PI*f0*f0*(k*T-t0)*(k*T-t0));
	}
	

	int s[X0][Z0]; //震源位置
	for(i=0;i<X0;i++)
	{
		for(j=0;j<Z0;j++)
		{
			if(i==X0/2 && j==Z0/2)
				s[i][j]=1;
			else
				s[i][j]=0;
		}
	}
	
	
	int V[X0][Z0];  //地下速度结构
	for(i=0;i<X0;i++)
	{
		for(j=0;j<Z0;j++)
		{
			V[i][j]=3000;
		}
	}


	float p1[X0][Z0],p2[X0][Z0],p3[X0][Z0];   //定义三个时间片
	for(i=0;i<X0;i++)
	{
		for(j=0;j<Z0;j++)
		{
			p1[i][j]=0;
			p2[i][j]=0;
			p3[i][j]=0;
		}
	}
	

	for(k=0;k<N;k++)         //三个时间片推导
	{
		if(k==0)
			;
		else
		{
			for(i=0;i<X0;i++)          //时间片递推
				for(j=0;j<Z0;j++)
				{
					p1[i][j]=p2[i][j];
					p2[i][j]=p3[i][j];
				}

			for(i=XH/2;i<X0-XH/2;i++)    //差分
			{
				for(j=ZH/2;j<Z0-ZH/2;j++)
				{
					float a=p2[i+1][j]-2*p2[i][j]+p2[i-1][j];
					float b=p2[i][j+1]-2*p2[i][j]+p2[i][j-1];
					p3[i][j]=V[i][j]*V[i][j]*T*T*(a/(DX*DX)+b/(DZ*DZ))+2*p2[i][j]-p1[i][j]+s[i][j]*f[k];
				}
			}
		}

		if(k==1999)
		{
			FILE *fp;
			if((fp=fopen("zhongxin2000.dat","wb"))!=NULL)
			{
				for(i=0;i<X0;i++)
					for(j=0;j<Z0;j++)
						fwrite(&p3[i][j],sizeof(float),1,fp);
			}
			fclose(fp);
		}
	}

	return 0;
}
