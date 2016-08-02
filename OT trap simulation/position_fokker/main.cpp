#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include<fstream>
#include<limits>

using namespace std;

//public constants:
const double pi=3.141592654;
const double kb=1.38*1e-23;//SI unit
const double temp=298; //(kelvin)
const double L=50*1e-9;

//public variables:
int N;
double Tmsr,dts,fs;
double k,r;
double D1,D2;
double D,G,fc;
string fileName;
int ave_num;

double KTP;


void initial_val();
void Langevin_position(double ,double ); //(solve langevin equation) calculates positions in an optical trap in presence of the thermal noise
double rand_normal(double, double ); //generates random number with normal (gaussian) distribution

//gets initial values
void initial_val(){

    cout <<"enter measurement time(sec):";
    cin>>Tmsr;

    cout<<"enter sampling frequency(Hz):";
    cin>>fs;

    cout<<"enter trap stiffness(pN/um):";
    cin>>k;
    k=k*1e-6; //SI

    cout<<"enter bead radius(um):";
    cin>>r;
    r=r*1e-6; //SI

    //calculate constants:
    G=6*pi*1e-3*r; //all units SI
    fc=k/(2.0*pi*G);
    D=kb*temp/G;
    N=Tmsr*fs;
    dts=double(1.0)/fs;
    D1=(-1.0)*2.0*pi*fc;//position is included in the main formula
    D2=2.0*D;

    KTP=kb*temp/(50*1e-9);
    double KTk;
    KTk=kb*temp/k;


    double A;
    A=D/(2*pi*pi);
    cout<< "G:"<<G<<" ,"<<"fc:"<<fc<<" ,"<<"D:"<<D<<" ,"<<"N:"<<N<<" ,A:"<<A<<"  T/dts:"<<Tmsr/dts<<"  kt/k:"<<KTk<< "  "<<sqrt(KTk)<<endl;
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

void Langevin_position(double *time,double *position1)
{
    time[0]=0.0;
    position1[0]=0.0;

    for(int n=0;n<N-1;n++)
    {
        time[n+1]=time[n]+dts;
        position1[n+1]=position1[n]+D1*dts*position1[n]+(sqrt(D2*dts)*rand_normal(0.0, 1.0));
    }
}

double forceWLC(double t,double xt, double &extension,double &xp)
{
    double F;
    double Vpull=50*1e-9/10000; //  nm/sec
    xp=Vpull*t;
    extension=xp-xt;

    F=(KTP)*((1/(4.0*(1.0-(extension/L))*(1.0-(extension/L))))+(extension/L)-0.25);
    return F;
}


void Langevin_position_WLC(double timeWLC,double positionWLC, double force, double ex)
{
    timeWLC=0.0;
    positionWLC=0.0;
    double xp=0.0;
    ofstream ofile;
    ofile.open("/Users/Naanaa/Desktop/test_WLCnoise1.csv");

    do{
            timeWLC=timeWLC+dts;
            positionWLC=positionWLC+(forceWLC(timeWLC,positionWLC,ex,xp)*dts/G)+D1*dts*positionWLC+(sqrt(D2*dts)*rand_normal(0.0, 1.0));
//            extension[n]=ex;
            ofile<<timeWLC<<","<<positionWLC<<","<<ex<<","<<k*positionWLC<<","<<xp<<endl;
    }while(timeWLC<Tmsr);
    ofile.close();
}



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
          arg = 2.0 * pi * (double)(k) / (double)(N);
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


int main()
{
    srand (time(NULL));

    initial_val();

    ofstream rand;
    rand.open("/Users/Naanaa/Desktop/WLCpulling/fs10k,vp50.csv");

//    ofstream deflection;
//    deflection.open("/Users/Naanaa/Desktop/comparisions/def,fc3k.txt");

//    double position1[N],time[N];
//    Langevin_position(time,position1);

    double timeWLC,positionWLC,force,extension;
//    double timeWLC[N],positionWLC[N],force[N],extension[N];
    Langevin_position_WLC(timeWLC,positionWLC, force, extension);

//    for(int i=0;i<N;i++)
//    {
//        rand<<timeWLC[i]<<","<<positionWLC[i]<<","<<force[i]<<","<<extension[i]<<endl;
//    }

//    double frequency[N/2],PS[N/2];
//    dft(position1,frequency,PS);

//    double stdev1;
//    stdev1=0.0;

//    for(int i=0;i<N;i++)
//    {
//        if(i<N/2)
//        {
//            rand<<time[i]<<","<<position1[i]<<","<<frequency[i]<<","<<PS[i]<<endl;
//        }
//        else
//        {
//            rand<<time[i]<<","<<position1[i]<<endl;
//        }
//        stdev1=stdev1+(position1[i]*position1[i]);
//    }

//    stdev1=stdev1/N;

//    rand<<endl;
//    rand<<"standard deviation: "<< sqrt(stdev1)<<","<<"trap stiffness:"<<kb*temp/stdev1<<endl;

//    cout <<"k1="<<k*1e6<<","<<"k2="<<(kb*temp*1e6/stdev1)<<endl<<"variance:"<<stdev1<<",standard deviation:" <<sqrt(stdev1)<<endl;
//    deflection<< "k1="<<k*1e6<<endl<<"k2="<<(kb*temp*1e6/stdev1)<<endl<<"variance(nm^2):"<<stdev1*1e18<<endl<<"standard deviation(nm^2):" <<sqrt(stdev1)*1e9<<endl;

    rand.close();
 //   deflection.close();

    return 0;
}
