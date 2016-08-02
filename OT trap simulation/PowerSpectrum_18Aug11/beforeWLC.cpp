#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <sstream>

using namespace std;

//public constants:
const double pi=3.141592654;
const double kb=1.38*1e-23;//SI unit
const double temp=298; //(kelvin)

//public variables:
double start_time = static_cast<double>(clock() / CLOCKS_PER_SEC);
int N;
double Tmsr,fs,dts,fnyq,dt; //dt time step between computations
double k; //trap stiffness
double c,dx;
double D,G,fc;
double D1,D2;
int ave_num;
string fileName;

void initial_val();
void Langevin_position_aliased(double ,double ); // calculates positions from Aliased langevin formula
void Langevin_position_dicrete_Fokker(double,double);
double rand_normal(double, double ); //generates random number with normal (gaussian) distribution
void dft(double, double, double);
//void theoretical_PS(int ,double ,double,double); //Pk theory
void integrate(double, double);

//calculates time in seconds
string current_time() {
    stringstream ss;
    double time = static_cast<double>(clock() / (double)CLOCKS_PER_SEC);
    ss.precision(3);
    ss << "[" << time-start_time << "] " ;
    return ss.str();
}

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

//Langevin position, Aliased Lorentzian as PS
void Langevin_position_aliased(double *time, double *position)
{
    double t=0.0; //  <time[j]<=t<=time[j+1]
    dt=dts/100.0;
    cout<<"dt for integral:"<<dt<<endl;

    time[0]=0.0;
    position[0]=0.0;

    double etha[N];
    double coef_etha;
    coef_etha=sqrt((4*pi*fc)/(1-(c*c)));
    double coef_exp=(-2)*pi*fc;
    double integ;

    for(int i=0;i<N;i++)
    {
        etha[i]=0;
    }

    for(int j=0;j<N-1;j++)
    {
        time[j]=j*dts;
        time[j+1]=(j+1)*dts;
        t=time[j];
        for(int i=0;i<=100;i++)
        {
            integ=(rand_normal(0.0,1.0)*exp(coef_exp*(time[j+1]-t)));
            if(i==0||i==100)
            {
                integ*=2;
            }
            etha[j]+=integ;
            t+=dt;
        }
        position[j+1]=c*position[j]+dx*coef_etha*etha[j]*dt;
        if(j%100==99){
          //  cout << current_time() << j<< " "<< c*position[j]<<"  "<<integ<<","<<dx*coef_etha*etha[j]*dt<< endl;
        }
    }
}


//langevin_discretized (Foker Plank)
void Langevin_position_dicrete_Fokker(double *time,double *position)
{
    time[0]=0.0;
    position[0]=0.0;

    for(int n=0;n<N-1;n++)
    {
        time[n+1]=time[n]+dts;
        position[n+1]=position[n]+D1*dts*position[n]+(sqrt(D2*dts)*rand_normal(0.0, 2.0));
        if(n%100==99)
        cout<<position[n]<<endl;
    }
}


//discrete Fourier Transform
void dft(double position_in[], double *f, double *mag2)
{
    long int k, j;
    double arg;
    double cosarg,sinarg;
    double out_real[N/2];
    double out_imag[N/2];

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
          out_imag[k-1]=out_imag[k-1]*dts;
          out_real[k-1]=out_real[k-1]*dts;
          f[k-1]=k*(1/Tmsr);
          mag2[k-1]=((out_real[k-1]*out_real[k-1])+(out_imag[k-1]*out_imag[k-1]))/Tmsr;// ith component of power spectrum
      }
}

//generates theoretical power spectra (Aliased lorentzian and lorentzian)
void theoretical_PS(int k,double &freq,double &Pk_Aliased, double &Pk_Lorentzian) //Pk theory
{
        freq=k/Tmsr;
        Pk_Aliased=(dx*dx)*dts/(1+(c*c)-(2*c*cos(2*pi*k/N))); //discretized (finite sampling) langevin
        Pk_Lorentzian=D/(2*pi*pi*((fc*fc)+(freq*freq))); //lorentzian
}

int main()
{
    srand (time(NULL));

    ofstream outPS;
    outPS.open("/Users/Naanaa/Desktop/PS_Data/fs100k_trap100.csv");

    initial_val(); //get initial values and calculates constants

    double position_Aliased[N]; //position from Aliased formula #15 in the Sorensen paper
    double position_Fokker[N];
    double time_pAliased[N]; //time from Aliased position (p) calculations
    double time_Fokker[N];
    double PS_Aliased[N/2];
    double Pk_Aliased[N/2];
    double PS_Fokker[N/2];
    double Pk_Fokker[N/2];
    double freq_al[N/2];
    double freq_fo[N/2];

    for(int i=0;i<N/2;i++)
    {
        PS_Fokker[i]=PS_Aliased[i]=0;
    }


    for(int i=0;i<ave_num;i++)
    {

//        Langevin_position_aliased( time_pAliased, position_Aliased);
//        dft(position_Aliased,freq_al,Pk_Aliased);

        Langevin_position_dicrete_Fokker(time_Fokker,position_Fokker);
        dft(position_Fokker,freq_fo,Pk_Fokker);

        for(int k=0;k<N/2;k++)
        {
//            PS_Aliased[k]=PS_Aliased[k]+Pk_Aliased[k];
            PS_Fokker[k]=PS_Fokker[k]+Pk_Fokker[k];

        if(i==ave_num-1)
        {
 //           PS_Aliased[k]/=ave_num;
            PS_Fokker[k]/=ave_num;

        }
    }
    }


    double freq_theory;
    double Pk_Aliased_theory,Pk_Lorentzian_theory;


    for(int k=1;k<=N;k++)
    {
        if(k<N/2)
        {
            //        outPS<<","<<time_pAliased[k-1]<<","<<position_Aliased[k-1]<<","<<freq_al[k-1]<<","<<Pk_Aliased[k-1]<<",";
                    outPS<<","<<time_Fokker[k-1]<<","<<position_Fokker[k-1]<<","<<freq_fo[k-1]<<","<<PS_Fokker[k-1]<<",";

                    theoretical_PS(k,freq_theory,Pk_Aliased_theory,Pk_Lorentzian_theory);
                    outPS<<","<<freq_theory<<","<<Pk_Aliased_theory<<","<<Pk_Lorentzian_theory<<endl;
        }
        else{
            outPS<<","<<time_Fokker[k-1]<<","<<position_Fokker[k-1]<<endl;

        }
    }

    outPS.close();

    return 0;
}

