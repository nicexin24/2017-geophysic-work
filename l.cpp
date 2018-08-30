
#include<iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
#define PI 3.1415
//#define N 2000     //时间片数量
#define T 0.001    //每个时间片的时间间隔
//#define t0 0   //时间延迟
int main()
{
//定义RICHER子波
cout<<"输入时间片："<<endl;
int N,i;
cin>>N;
cout<<"输入子波频率个数："<<endl;
int h;
cin>>h;
cout<<"输入子波频率："<<endl;
int ff[h];
for(i=0;i<h;i++)
cin>>ff[i];
cout<<"时间延迟"<<endl;
float t0;
cin>>t0;
	float f[h][N]; 
int k;
	int A=1;    //振幅
	for(i=0;i<h;i++)
	for(k=0;k<N;k++)
	{
		f[i][k]=A*exp(-PI*PI*ff[i]*ff[i]*pow(k*T-t0,2))*(1-2*PI*PI*ff[i]*ff[i]*pow(k*T-t0,2));

	}
                   FILE *fp;
			if((fp=fopen("l.dat","wb"))!=NULL)
			{
				for(int i=0;i<h;i++)
                                   for(k=0;k<N;k++)
						fwrite(&f[i][k],sizeof(float),1,fp);
			}
                         else 
                      cout<<"fail to build new file";
			fclose(fp);
}
