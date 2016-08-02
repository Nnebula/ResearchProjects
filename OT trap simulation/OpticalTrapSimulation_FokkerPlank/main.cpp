#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include<fstream>
#include <time.h>

using namespace std;

int N;
double T,dt;
double k;
double D1,D2;
double D,G,fc;

const double pi=3.141592654;
const double kb=1.38*1e-23;//SI unit
const double temp=298; //(kelvin)


void initial_val();
void position_cal(double ,double );
double rand_normal(double, double );
void Langevin_WLC(double,double);

//gets initial values
void initial_val()
{
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

    //calculate constants:
    G=6*pi*1e-3*r; //all units SI
    fc=k/(2.0*pi*G);
    D=kb*temp/G;
    N=T/dt;

    D1=(-1.0)*2.0*pi*fc;//position is included in the main formula
    D2=2.0*D;
    cout<< "G:"<<G<<" ,"<<"fc:"<<fc<<" ,"<<"D:"<<D<<" ,"<<"N:"<<N<<endl;
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

//solve langevin equation for a trapped bead
void position_cal(double *time,double *position)
{
    ofstream LangevinFile;
    LangevinFile.open("/Users/Naanaa/Desktop/randomNumber/position_rand_compare.csv");

    time[0]=0.0;
    position[0]=0.0;

    for(int n=0;n<N-1;n++)
    {
        LangevinFile<<time[n]<<","<<position[n]<<","<<rand_normal(0.0,2.0)<<","<<rand_normal(0.0,1.0)<<endl;
        position[n+1]=position[n]+D1*dt*position[n]+(sqrt(D2*dt)*rand_normal(0.0, 2.0));        
		time[n+1]=time[n]+dt;
    }
    LangevinFile.close();
}



////////////
int main()
{
    srand (time(NULL));

    initial_val(); //get initial values and calculates constants

    double position[N];
    double time[N];
    position[0]=0;
    time[0]=0.0;

//    position_cal(time,position);//Langevinequation solution (fokker_plank)
//    dft(position, FT_pos_real, FT_pos_imag, f, Pk); //Discrete fourier transform
//    WLC(extension, force);

//    for(int i=0;i<wlc;i++)
//    {
//        Langevin_WLC(time[i],position_WLC[i]);
//    }
    position_cal(time,position);

    return 0;
}
