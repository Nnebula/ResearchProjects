#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include<fstream>
#include <time.h>

using namespace std;

//public constants:
const double pi=3.141592654;
const double kb=1.38*1e-23;//SI unit
const double temp=298; //(kelvin)

//public variables:
int N;
double T,dt;
double k;
double D1,D2;
double D,G,fc;
string fileName;
int ave_num;

void initial_val();
void Langevin_position(double ,double ); //(solve langevin equation) calculates positions in an optical trap in presence of the thermal noise
double rand_normal(double, double ); //generates random number with normal (gaussian) distribution
void dft(double, double, double);

//gets initial values
void initial_val(){

    cout <<"enter measurement time(sec):";
    cin>>T;

    cout<<"enter dt(sec):";
    cin>>dt;

    cout<<"enter trap stiffness(pN/um):";
    cin>>k;
    k=k*1e-6; //SI

    double r;
    cout<<"enter bead radius(um):";
    cin>>r;
    r=r*1e-6; //SI

    cout<<"enter file name: date/number:";
    cin>>fileName;

    cout<<"enter the number of Power Spectra to be averaged on:";
    cin>>ave_num;

    //calculate constants:
    G=6*pi*1e-3*r; //all units SI
    fc=k/(2.0*pi*G);
    D=kb*temp/G;
    N=T/dt;

    D1=(-1.0)*2.0*pi*fc;//position is included in the main formula
    D2=2.0*D;

    double A;
    A=D/(2*pi*pi);
    cout<< "G:"<<G<<" ,"<<"fc:"<<fc<<" ,"<<"D:"<<D<<" ,"<<"N:"<<N<<" ,A:"<<A<<endl;
}

//generates gaussian random number
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

//solves langevin equation (no other force than the trap and thermal noise)
void Langevin_position(double *time,double *position)
{
    ofstream LangevinFile;
    LangevinFile.open(("/Users/Naanaa/Documents/Programming/Research/OT trap simulation/PowerSpectrum/data/LangevinPosition_" +fileName+ ".csv").c_str());

    time[0]=0.0;

    for(int n=0;n<N-1;n++)
    {
        time[n+1]=time[n]+dt;
        position[n+1]=position[n]+D1*dt*position[n]+(sqrt(D2*dt)*rand_normal(0.0, 2.0));
    }

//    LangevinFile<<time[n]<<","<<position[i][n]<<endl;
    LangevinFile.close();
}

//discrete Fourier Transform
void dft(double position_in[], double *k, double *mag2)
{
  long int i, j;
  double arg;
  double cosarg,sinarg;
  double out_real[N];
  double out_imag[N];

  for(i=0; i<N; i++)
      {
        out_real[i] = 0.0;
        out_imag[i] = 0.0;
        arg = 2.0 * pi * (double)i / (double)N;
        for(j=0; j<N; j+=1)
        {
          cosarg = cos(j * arg);
          sinarg = sin(j * arg);
          out_imag[i] += (position_in[j] * sinarg);
          out_real[i] += (position_in[j] * cosarg);
        }
        out_imag[i]=out_imag[i]*dt;
        out_real[i]=out_real[i]*dt;
        k[i]=i*(1/T);
        mag2[i]=((out_real[i]*out_real[i])+(out_imag[i]*out_imag[i]))/T;// ith component of power spectrum
    }
}

int main()
{
    srand (time(NULL));

    initial_val(); //get initial values and calculates constants

    double position[N];
    double time[N];
    time[0]=0.0;
    double Pk[N];
    double f[N];
    double PS[N];
    ofstream outPS;
    outPS.open(("/Users/Naanaa/Documents/Programming/Research/OT trap simulation/PowerSpectrum/data/PowerSpec_" +fileName+ ".csv").c_str());

    for(int i=0;i<N;i++)
    {
        PS[i]=0.0;
    }

    for(int i=0;i<ave_num;i++)
    {
        Langevin_position(time,position);//Langevinequation solution (fokker_plank)
        dft(position, f, Pk); //Discrete fourier transform

        for(int k=0;k<N;k++)
        {
            PS[k]=Pk[k]+PS[k];
            if(i==ave_num-1)
            {
                PS[k]=PS[k]/ave_num;
                outPS<<f[k]<<","<<PS[k]<<endl; //average power spectrum over 100 PS
            }
        }
    }
    outPS.close();

    return 0;
}



