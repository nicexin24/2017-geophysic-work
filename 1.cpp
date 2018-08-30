
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
#define PI 3.1415
#define N 2000     //时间片数量
#define T 0.001    //每个时间片的时间间隔
//#define t0 0.05   //时间延迟
#define X 350      //x单元数
#define Z 350      //z单元数
#define DX 10    //x轴步长
#define DZ 10  //z轴步长
#define M  5     //2M为空间差分阶数
//定义数组长度
const int X0=X+2*M;     
const int Z0=Z+2*M; 

int i,j,k;
//差分系数 
float cc[5][5]={1.0000000e+00,0.0000000e+00,0.0000000e+00,0.0000000e+00,0.0000000e+00,
   1.3333333e+00,-8.3333333e-02,0.0000000e+00,0.0000000e+00,0.0000000e+00,
   1.5000000e+00,-1.5000000e-01,1.1111111e-02,0.0000000e+00,0.0000000e+00,
   1.6000000e+00,  -2.0000000e-01,2.5396825e-02,-1.7857143e-03,0.0000000e+00,
   1.6666667e+00,-2.3809524e-01,3.9682540e-02,-4.9603175e-03,3.1746032e-04};
//定义时间片
	  float p1[X0][Z0],p2[X0][Z0],p3[X0][Z0];   

int main()
{
//输入参数
cout<<"输入时间节点："<<endl;
int k0;
cin>>k0;
cout<<"输入子波频率："<<endl;
int f0;
cin>>f0;
cout<<"时间延迟"<<endl;
float t0;
cin>>t0;
cout<<"输入类型vsp（1）or普通反射(0):"<<endl;
int r;
cin>>r;
/*cout<<"输入查分阶数："<<endl;
int M;
cin>>M;*/
float u[k0][Z0];
//计算具体差分系数
float c[M];
i=M-1;
for (j=0;j<M;j++)
c[j]=cc[i][j];
       
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
			if(i==X0/2 && j==M)
				s[i][j]=1;
			else
				s[i][j]=0;
		}
	//地下速度结构
cout<<"选择速度结构："<<endl<<"均匀介质 0"<<endl<<"二层介质 1"<<endl<<"三层介质 2"<<endl<<"斜层介质 3"<<endl<<"阶梯介质 4"<<endl<<"斜层&多层介质 5"<<endl<<"上低速&下斜层介质 6"<<endl<<"高速球体介质 7"<<endl;
int h;
cin>>h;
	int V[X0][Z0];  
	for(i=0;i<X0;i++)
		for(j=0;j<Z0;j++)
		{
switch (h){
case 0:V[i][j]=3000; break;                   
case 1:{if (j>150)V[i][j]=5000; else V[i][j]=3000;}break;
case 2:{if (j<100) V[i][j]=3000;else  if(j<150) V[i][j]=4000;else V[i][j]=5000;}break;
case 3:{if(j>0.8*i) V[i][j]=5000;else  V[i][j]=3000;}break;
case 4:{if((j>100&&i<X0/2) ||(j>150&&i>=X0/2))V[i][j]=5000;else V[i][j]=3000; }break;
case 5:{if(j>0.4*i&&j<200) V[i][j]=4000;else if(j>=200) V[i][j]=5000;else V[i][j]=3000; }break;
case 6:{if(j<100) V[i][j]=3000;else if(j>Z0-1.5*i) V[i][j]=5000;else V[i][j]=4000;}break;
case 7:{float vv=pow(i-X0/2,2)+pow(j-100,2);if (vv<2500)V[i][j]=5000; else V[i][j]=3000;}break;
		}
}
		
		
      
//赋初值
	for(i=0;i<X0;i++)
		for(j=0;j<Z0;j++)
		{
			p1[i][j]=0;p2[i][j]=0;p3[i][j]=0;
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
				}
		//划分空间区域
		
				//正常区域
				
				for(i=M;i<X+M+1;i++)          
				      for(j=M;j<Z+M+1;j++) 
				{
					float a(0),b(0);
					//计算空间导数算子
					for (int mm=1;mm<M+1;mm++)
					{int m=mm;a=a+c[m-1]*(p2[i-m][j]-2*p2[i][j]+p2[i+m][j]);
						b=b+c[m-1]*(p2[i][j-m]-2*p2[i][j]+p2[i][j+m]);
					}	
	                   p3[i][j]=pow(V[i][j],2)*T*T*(a/(DX*DX)+b/(DZ*DZ))+2*p2[i][j]-p1[i][j]+s[i][j]*f[k];
				}
                
		for (i=0;i<Z0;i++)
                       u[k-1][i]=p3[X0/2][i];
		
	
              if(k==k0)
			
		{
			FILE *fp;
			if((fp=fopen("1.dat","wb"))!=NULL)
			{
				if(r==1)
                                   {for(i=0;i<k0;i++)
					for(j=0;j<Z0;j++)
						fwrite(&u[i][j],sizeof(float),1,fp);}
else
                                    { for(i=0;i<X0;i++)
					for(j=0;j<Z0;j++)
						fwrite(&p3[i][j],sizeof(float),1,fp);}
			}
                         else 
                      cout<<"fail to build new file";
			fclose(fp);
                break;
		}
	}
	return 0;
}
