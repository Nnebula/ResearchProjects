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
double Tmsr,fs,dts,fnyq,dt; //dt time step between computations
double k; //trap stiffness
double c,dx;
double D,G,fc;
int ave_num;
string fileName;

void initial_val();
void Langevin_position_aliased(double ,double ); // calculates positions from Aliased langevin formula
double rand_normal(double, double ); //generates random number with normal (gaussian) distribution
void dft(double, double, double);
void theoretical_PS_notLorentzian(int ,double ,double,double); //Pk theory

//gets initial values
void initial_val(){

    cout <<"enter measurement time(sec) Tmsr:";
    cin>>Tmsr;

    cout<<"enter sampling frequency(Hz) fs:";
    cin>>fs;

    cout<<"enter trap stiffness(pN/um) k:";
    cin>>k;
    k=k*1e-6; //SI

    double r;
    cout<<"enter bead radius(um) r:";
    cin>>r;
    r=r*1e-6; //SI

    cout<<"enter file name: date/format:";
    cin>>fileName;

    cout<<"enter the number of Power Spectra to be averaged on:"; //ave_num
    cin>>ave_num;

    //calculate constants:
    G=6*pi*1e-3*r; //all units SI
    fc=k/(2.0*pi*G);
    D=kb*temp/G;
    N=Tmsr*fs; //number of data generated

    dts=1.0/fs;
    fnyq=fs/2.0;

    c=exp((-1)*pi*fc/fnyq);
    dx=sqrt((1-(c*c))*D/(2*pi*fc));

    D1=(-1.0)*2.0*pi*fc;//position is included in the main formula
    D2=2.0*D;

    double A;
    A=D/(2*pi*pi);
    cout<< "G:"<<G<<" , fc:"<<fc<<" , D:"<<D<<" , A:"<<A<<endl;
    cout<<"fs:"<<fs<<" , dts:"<<dts<<" , N:"<<N<<" , c:"<<c<<" , dx:"<<dx<<endl;
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

//NonAliased Position
void position_cal(double *time,double *position)
{
    ofstream LangevinFile;
    LangevinFile.open("/Users/Naanaa/Documents/Programming/Research/OT trap simulation/OpticalTrapSimulation_FokkerPlank/data/Langevin_WLC.csv");

    time[0]=0.0;
    position[0]=0.0;

    for(int n=0;n<N-1;n++)
    {
        LangevinFile<<time[n]<<","<<position[n]<<endl;
        time[n+1]=time[n]+dt;
        position[n+1]=position[n]+D1*dt*position[n]+dt*(WLC_force(position[n])/G)+(sqrt(D2*dt)*rand_normal(0.0, 2.0));
    }
    LangevinFile.close();
}

//aliased lorentzian equation (position)
void Langevin_position_aliased(double *time, double *position)
{
    double t=0.0; //  <time[j]<=t<=time[j+1]
    dt=dts/10.0;
    cout<<"dt for integral:"<<dt<<endl;

    time[0]=0.0;
    position[0]=0.0;

    ofstream LangevinFile;
    LangevinFile.open(("/Users/Naanaa/Documents/Programming/PowerSpectrum_18Aug11/data/Position_Aliased_" +fileName+ ".csv").c_str());

    double etha[N];
    for(int i=0;i<N;i++)
    {
        etha[i]=0;
    }

    double coef_etha;
    coef_etha=sqrt((4*pi*fc)/(1-(c*c)));

    for(int j=0;j<N-1;j++)
    {
        time[j]=j*dts;
        time[j+1]=(j+1)*dts;
        t=time[j];
        do
        {
            etha[j]=etha[j]+(rand_normal(0,1)*dt*exp((-2)*pi*fc*(time[j+1]-t)));
            if(j%100==99)
                cout<<time[j] <<"  "<<t<<endl;
            t=t+dt;
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

  ofstream outPS;
  outPS.open(("/Users/Naanaa/Documents/Programming/PowerSpectrum_18Aug11/data/Position_Aliased_" +fileName+ ".csv").c_str());


  for(k=1; k<=N/2; k++)
      {
        out_real[k-1] = 0.0;
        out_imag[k-1] = 0.0;
        arg = 2.0 * pi * (double)(k) / (double)N;
        for(j=0; j<N; j+=1)
        {
          cosarg = cos(j * arg);
          sinarg = sin(j * arg);
          out_imag[k-1] += (position_in[j] * sinarg);
          out_real[k-1] += (position_in[j] * cosarg);
        }
        out_imag[k-1]=out_imag[k-1]*dt;
        out_real[k-1]=out_real[k-1]*dt;
        f[k-1]=k*(1/Tmsr);
        mag2[k-1]=((out_real[k-1]*out_real[k-1])+(out_imag[k-1]*out_imag[k-1]))/Tmsr;// ith component of power spectrum
        outPS<<f[k-1]<<","<<mag2[k-1]<<endl;
    }
  outPS.close();
}

//generates power spectrum of non-lorentzian form(due to finite sampling) and lorentzian
void theoretical_PS_notLorentzian(int k,double &freq,double *Pk_t, double *Pk_Lor) //Pk theory
{
        freq=k/Tmsr;
        Pk_t[k]=(dx*dx)*dts/(1+(c*c)-(2*c*cos(2*pi*k/N))); //discretized (finite sampling) langevin
        Pk_Lor[k]=D/(2*pi*pi*((fc*fc)+(freq*freq))); //lorentzian
}

int main()
{
    srand (time(NULL));

    initial_val(); //get initial values and calculates constants

    double position_Aliased[N]; //position from Aliased formula #15 in the Sorensen paper
    double time_pAliased[N]; //time from Aliased position (p) calculations
    double freq[N/2];
    double PS_Aliased[N/2];

    Langevin_position_aliased( time_pAliased, position_Aliased);
    dft(position_Aliased,freq,PS_Aliased);

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
                outPS<<f[k]<<","<<PS[k]<<endl; //average power spectrum over ave_num Power spectra
            }
        }
    }

//    double fr;
//    double Pkt[N/2]; //non-lorentzian
//    double Pk_Lor[N/2]; // lorentzian
//    for( int k=1;k<=N/2;k++)
//    {
//        theoretical_PS_notLorentzian(k,fr,Pkt,Pk_Lor);
//        outPS<<k<<","<<fr<<","<<Pkt[k]<<","<<Pk_Lor[k]<<endl;
//    }

   //


    return 0;
}




