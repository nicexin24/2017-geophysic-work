#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
#define h 10    //定义空间坐标的标准步长，两个方向采取相同步长
#define X 500    //定义空间x方向所取点数
#define Z 500    //定义空间Z方向所取点数
#define XW 10   //定义空间x方向延伸长度
#define ZW 10    //定义空间Z方向延伸长度
#define T 2     //定义研究的时间长度
#define N 2000  //时间片个数
#define ti 0.001  //定义时间间隔
#define pi 3.1415
int main()
{
 float p1[X][Z],p2[X][Z],p3[X][Z],t;
 int i,j,k;
 //定义RICHER子波
 float f[N];    
	int f0=30;  //主频
	int A=1;    //振幅
 for(k=0;k<N;k++)
	{
		f[k]=A*exp(-pi*pi*f0*f0*k*T*k*T)*(1-2*pi*pi*f0*f0*k*T*k*T);
	}


  //地下速度赋值
 int V[X][Z];       
	for(i=0;i<X;i++)
	{
		for(k=0;k<Z;k++)
		{
			V[i][k]=3000;
		}
	}

	 //震源位置空间函数
	int s[X][Z];
	for(i=0;i<X;i++)
	{
		for(k=0;k<Z;k++)
		{
			if(i==X/2 && k==Z/2)
				s[i][k]=1;
			else
				s[i][k]=0;
		}
	}
	//p1，p2，p3赋初值
	for (i=0;i<500;i++)
	 for (k=0;k<500;k++)
		{ p1[i][k]=0;
         p2[i][k]=0;
		 p3[i][k]=0;
	     }
	//p3时间片震动情况的推导
	//利用波动方程，差分计算
	
 for (j=0;j<2000;j++)
 {
	 //三个时间片的循环递推
	 for (i=0;i<500;i++)
	 for (k=0;k<500;k++)
	 {p1[i][k]=p2[i][k];p2[i][k]=p3[i][k];}

	 //计算新的p3
	 for(i=XW/2;i<X-XW/2;i++)    //差分
				for(k=ZW/2;j<Z-ZW/2;k++)
				{
					float a=p2[i+1][k]-2*p2[i][k]+p2[i-1][k];
					float b=p2[i][k+1]-2*p2[i][k]+p2[i][k-1];
					p3[i][k]=V[i][k]*V[i][k]*T*T*(a/(h*h)+b/(h*h))+2*p2[i][k]-p1[i][k]+s[i][k]*f[j];
				}	
      if(j==1999)
		{
			FILE *fp;
			if(fopen_s(&fp,"0522.dat","wb")==0)	
			{
				for(i=0;i<X;i++)
					for(k=0;k<Z;k++)
						fwrite(&p3[i][k],sizeof(float),1,fp);
			}
			fclose(fp);
		}
 }

return 0;
}