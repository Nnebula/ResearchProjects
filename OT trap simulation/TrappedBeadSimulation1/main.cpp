
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include<fstream>
#include <time.h>

using namespace std;

double k, T, G, dt,t;
const long int length=10000;
const double pi=3.14;

void get_k();
void get_T();
void get_G();
void get_dt();

double rand_normal(double, double);
void calculate(double);
void PS_theory(double);
void PS_exp(double,double,double);


//get initial values:

void get_k()
{
    cout << "enter trap stiffness(pN/um):";
    cin>> k;
}


void get_T()
{
    cout<<"enter the temprature(centigrade):";
    cin >> T;
    T=(T+273)*1.38*1e-5;//units: pN.um
}


void get_G()
{
    cout<< "enter the trapped bead radius in um:";
    cin >> G;
    G=6*3.14*G*1e+3; //units: um^2/sec
}

void get_dt()
{
    cout<< "enter the desired dt for calculations in seconds:";
    cin>>dt;
}


//gaussian random number generator (box-muller transfor),
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


// calculate position of the bead in the trap at each time

void calculate(double *position)
{
    double t=0.0;

    for(int i=0;i<length;i++)
    {
        position[i]=0;
    }

    for(int i=0;i<length;i++)
    {
        position[i]=position[i]+(-(k/G)*position[i]+sqrt(2*T/G)*rand_normal(0.0,1.0))*dt;
       t=t+dt;
   }
}



//Discrete Fourier Transform
void dft(long int length, double input_sample[], double *out_imag, double *out_real)
{
  long int i, j;
  double arg;
  double cosarg,sinarg;

  for(i=0; i<length; i++)
  {
    out_real[i] = 0;
    out_imag[i] = 0;
    arg = -1.0 * 2.0 * 3.141592654 * (double)i / (double)length;
    for(j=0; j<length; j+=1)
    {
      cosarg = cos(j * arg);
      sinarg = sin(j * arg);
      out_imag[i] += (input_sample[j] * sinarg);
      out_real[i] += (input_sample[j] * cosarg);
    }
    out_imag[i]=out_imag[i]*dt;
    out_real[i]=out_real[i]*dt;
  }
}

//power spectrum calculation (experimental)
void PS_exp(double out_real[],double out_imag[],double *Pk)
{
    double mag2[length];
    for(int k=0;k<length;k++)
    {
        mag2[k]=(out_real[k]*out_real[k])+(out_imag[k]*out_imag[k]);
        Pk[k]=mag2[k]/length*dt;
    }
}

//power spectrum calculation (theory)

void PS_theory(double *out_PS)
{
    double fc=k/(2*3.14*G);
    double fk;
    double D=T/G;

    for(int i=0;i<length;i++)
    {
        out_PS[i]=0;
    }

    for(int k=0;k<length;k++)
    {
        fk=k/dt;
        out_PS[k]=(D/(2*pi*pi))/((fk*fk)+(fc*fc));
    }
}


int main()
{
     ofstream FTfile;
     srand(time(NULL));
     get_k();
     get_T();
     get_G();
     get_dt();

     double out_real[length],out_imag[length];
     double position[length];
     double out_PS[length];
     double Pk[length];

     calculate(position);
     dft(length, position,out_real,out_imag);
     PS_exp(out_real,out_imag,Pk);
     PS_theory(out_PS);

     FTfile.open("/Users/Naanaa/Desktop/OT_simulations/FT.csv");

     double t=0;

     for(int j=0;j<length;j++)
     {
         FTfile<<t<<","<<position[j]<<","<<out_real[j]<<","<<out_imag[j]<<","<<out_PS[j]<<","<<","<<j<<","<<Pk[j]<<endl;
         t=t+dt;
     }

     FTfile.close();

    return 0;
}
