#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;
#define h 10    //����ռ�����ı�׼���������������ȡ��ͬ����
#define X 500    //����ռ�x������ȡ����
#define Z 500    //����ռ�Z������ȡ����
#define XW 10   //����ռ�x�������쳤��
#define ZW 10    //����ռ�Z�������쳤��
#define T 2     //�����о���ʱ�䳤��
#define N 2000  //ʱ��Ƭ����
#define ti 0.001  //����ʱ����
#define pi 3.1415
int main()
{
 float p1[X][Z],p2[X][Z],p3[X][Z],t;
 int i,j,k;
 //����RICHER�Ӳ�
 float f[N];    
	int f0=30;  //��Ƶ
	int A=1;    //���
 for(k=0;k<N;k++)
	{
		f[k]=A*exp(-pi*pi*f0*f0*k*T*k*T)*(1-2*pi*pi*f0*f0*k*T*k*T);
	}


  //�����ٶȸ�ֵ
 int V[X][Z];       
	for(i=0;i<X;i++)
	{
		for(k=0;k<Z;k++)
		{
			V[i][k]=3000;
		}
	}

	 //��Դλ�ÿռ亯��
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
	//p1��p2��p3����ֵ
	for (i=0;i<500;i++)
	 for (k=0;k<500;k++)
		{ p1[i][k]=0;
         p2[i][k]=0;
		 p3[i][k]=0;
	     }
	//p3ʱ��Ƭ��������Ƶ�
	//���ò������̣���ּ���
	
 for (j=0;j<2000;j++)
 {
	 //����ʱ��Ƭ��ѭ������
	 for (i=0;i<500;i++)
	 for (k=0;k<500;k++)
	 {p1[i][k]=p2[i][k];p2[i][k]=p3[i][k];}

	 //�����µ�p3
	 for(i=XW/2;i<X-XW/2;i++)    //���
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