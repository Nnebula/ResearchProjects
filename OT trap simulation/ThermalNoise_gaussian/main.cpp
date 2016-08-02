#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include<fstream>
#include <time.h>

using namespace std;


const long int length=100;
const int ave_num=100;
const double pi=3.14;

double rand_normal(double, double);
void dft_thermal(double, double, double);

//gaussian random number generator (box-muller transform),
 /*1 Copyright (c) 2011 the authors listed at the following URL, and/or
2 the authors of referenced articles or incorporated external code:
3 http://en.literateprograms.org/Box-Muller_transform_(C)?action=history&offset=20060711193108
4
5 Permission is hereby granted, free of charge, to any person obtaining
6 a copy of this software and associated documentation files (the
7 "Software"), to deal in the Software without restriction, including
8 without limitation the rights to use, copy, modify, merge, publish,
9 distribute, sublicense, and/or sell copies of the Software, and to
10 permit persons to whom the Software is furnished to do so, subject to
11 the following conditions:
12
13 The above copyright notice and this permission notice shall be
14 included in all copies or substantial portions of the Software.
15
16 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
17 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
18 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
19 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
20 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
21 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
22 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
23
24 Retrieved from: http://en.literateprograms.org/Box-Muller_transform_(C)?oldid=7011
25 */

double rand_normal(double mean, double stddev)
{
    static double n2 = 0.0;
    static int n2_cached = 0;

    if (!n2_cached)
    {

        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        } while (r == 0.0 || r > 1.0);

        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;

            double result = n1*stddev + mean;

            n2_cached = 1;
            return result;
        }

     }

    else
    {
            n2_cached = 0;
            return n2*stddev + mean;
    }
 }


void dft_thermal(double position_in, double *mag2, double *mag4)
{
    double arg;
    double cosarg,sinarg;
    double out_real[length];
    double out_imag[length];
    double f[length];


    for(int k=1; k<=length; k++)
        {
          out_real[k-1] = 0.0;
          out_imag[k-1] = 0.0;
          arg = 2.0 * pi * (double)(k) / (double)length;

          for(int j=0; j<length; j+=1)
          {
            cosarg = cos(j * arg);
            sinarg = sin(j * arg);
            out_imag[k-1] += (position_in * sinarg);
            out_real[k-1] += (position_in * cosarg);
          }
          out_imag[k-1]=out_imag[k-1]/length;
          out_real[k-1]=out_real[k-1]/length;
          f[k-1]=k*(1/length);
          mag2[k-1]=pow(out_real[k-1],2)+pow(out_imag[k-1],2);// ith component of power spectrum
          mag4[k-1]=pow(mag2[k-1],2);
      }
}



int main()
{
    double mean;
    double stddev;

    cout<<"enter mean:";
    cin>>mean;
    cout<<"enter width:";
    cin>>stddev;

    double in_thermal[ave_num];
    double FT_thermal_2[length],FT_thermal_4[length];
    double sum2[length],sum4[length];

    for(int i=0;i<length;i++)
    {
        sum2[i]=sum4[i]=0.0;
    }


    ofstream FTfile;
    srand(time(NULL));
    double x;

    for(int i=0;i<ave_num;i+=1)
    {
        in_thermal[i]=rand_normal(mean,stddev);
        in_thermal[i]=x;

        dft_thermal(x, FT_thermal_2,FT_thermal_4);
        for(int j=0;j<length;j++)
        {
            sum2[j]=FT_thermal_2[j]+sum2[j];
            sum4[j]=FT_thermal_4[j]+sum4[j];
        }
    }

    FTfile.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/Powerspectrum/thermalNoise_FT.txt");
    for(int k=0;k<ave_num;k+=1)
    {
        FTfile<< k <<","<<sum2[k]/ave_num<<","<<sum4[k]/ave_num<<endl;
    }
    FTfile.close();
    return 0;
}
