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
double T,fs,dts;
double k;
double c,dt,dx;
double D,G,fc;
string fileName;
int ave_num;

void initial_val();
void Langevin_position(double ,double ); //(solve langevin equation) calculates positions in an optical trap in presence of the thermal noise
double rand_normal(double, double ); //generates random number with normal (gaussian) distribution
void dft(double, double, double);
void theoretical_PS_notLorentzian(int ,double ,double,double); //Pk theory

//gets initial values
void initial_val(){

    cout <<"enter measurement time(sec):";
    cin>>T;

    cout<<"enter sampling frequency(Hz):";
    cin>>fs;

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
    N=T*fs;

    dts=1/fs;
    double fnyq=0.5*fs;

    c=exp((-1)*pi*fc/fnyq);
    dx=sqrt((1-(c*c))*D/(2*pi*fc));

    double A;
    A=D/(2*pi*pi);
   // cout<< "G:"<<G<<" ,"<<"fc:"<<fc<<" ,"<<"D:"<<D<<" ,"<<"N:"<<N<<" ,A:"<<A<<endl;
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
    LangevinFile.open(("/Users/Naanaa/Documents/Programming/PowerSpectrum_antiAliased/data/position_discretized_" +fileName+ ".csv").c_str());

    double t; //  <time[j]<=t<=time[j+1]
    dt=dts/10.0;
    cout<<dt<< "hello honey :*";
    time[1]=0.0;
    position[1]=0.0;

    double etha[N];
    for(int i=1;i<=N;i++)
    {
        etha[i]=0;
    }

    double coef_etha;
    coef_etha=sqrt((4*pi*fc)/(1-(c*c)));

    for(int j=1;j<=N-1;j++)
    {
        time[j+1]=time[j]+dts;
        t=time[j];
        do
        {
            etha[j]=etha[j]+(rand_normal(0,1)*dt*exp((-2)*pi*fc*(time[j+1]-t)));
            t=t+dt;
               // cout<< time[j]<<"  "<<t<<endl;
        }while(t<=time[j+1]);

    position[j+1]=c*position[j]+dx*coef_etha*etha[j];
    LangevinFile<<time[j]<<","<<position[j]<<endl;
    }

    LangevinFile.close();
}

//discrete Fourier Transform
void dft(double position_in[], double *f, double *mag2)
{
  long int k, j;
  double arg;
  double cosarg,sinarg;
  double out_real[N];
  double out_imag[N];

  for(k=1; k<=N; k++)
      {
        out_real[k] = 0.0;
        out_imag[k] = 0.0;
        arg = 2.0 * pi * (double)(k-N/2) / (double)N;
        for(j=0; j<N; j+=1)
        {
          cosarg = cos(j * arg);
          sinarg = sin(j * arg);
          out_imag[k] += (position_in[j] * sinarg);
          out_real[k] += (position_in[j] * cosarg);
        }
        out_imag[k]=out_imag[k]*dt;
        out_real[k]=out_real[k]*dt;
        f[k]=k*(1/T);
        mag2[k]=((out_real[k]*out_real[k])+(out_imag[k]*out_imag[k]))/T;// ith component of power spectrum
    }
}

//generates power spectrum of non-lorentzian form(due to finite sampling) and lorentzian
void theoretical_PS_notLorentzian(int k,double &freq,double *Pk_t, double *Pk_Lor) //Pk theory
{
        freq=k/T;
        Pk_t[k]=(dx*dx)*dts/(1+(c*c)-(2*c*cos(2*pi*k/N))); //discretized (finite sampling) langevin
        Pk_Lor[k]=D/(2*pi*pi*((fc*fc)+(freq*freq))); //lorentzian
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
    outPS.open(("/Users/Naanaa/Documents/Programming/PowerSpectrum_antiAliased/data/discretized_PS_" +fileName+ ".csv").c_str());

    for(int i=0;i<N;i++)
    {
        PS[i]=0.0;
    }

//    for(int i=0;i<ave_num;i++)
//    {
//        Langevin_position(time,position);//Langevinequation solution (fokker_plank)
//        dft(position, f, Pk); //Discrete fourier transform

//        for(int k=0;k<N;k++)
//        {
//            PS[k]=Pk[k]+PS[k];
//            if(i==ave_num-1)
//            {
//                PS[k]=PS[k]/ave_num;
//                outPS<<f[k]<<","<<PS[k]<<endl; //average power spectrum over ave_num Power spectra
//            }
//        }
//    }

    double fr;
    double Pkt[N/2]; //non-lorentzian
    double Pk_Lor[N/2]; // lorentzian
    for( int k=1;k<=N/2;k++)
    {
        theoretical_PS_notLorentzian(k,fr,Pkt,Pk_Lor);
        outPS<<k<<","<<fr<<","<<Pkt[k]<<","<<Pk_Lor[k]<<endl;
    }

    outPS.close();


    return 0;
}



