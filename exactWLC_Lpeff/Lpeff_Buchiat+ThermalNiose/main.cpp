
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include<fstream>
#include<limits>

using namespace std;

//public constants:
const double PI=3.141592654;
const double KB=1.38*1e-23;//SI unit
const double TEMP=298; //(kelvin)

//public variables:
double a2,a3,a4,a5,a6,a7;
double L;//contour length
int N;
double Tmsr,dts,fs;
double k,r;
double D1,D2;
double D,G,fc;

double KTk=KB*TEMP/k;
double KTP=KB*TEMP/(50*1e-9); //SI units
double KT=KB*TEMP; //SI units



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
    G=6*PI*1e-3*r; //all units SI
    fc=k/(2.0*PI*G);
    D=KB*TEMP/G;
    N=Tmsr*fs;
    dts=double(1.0)/fs;
    D1=(-1.0)*2.0*PI*fc;//position is included in the main formula
    D2=2.0*D;
    L=3.958*1e-6; //SI units


    a2=-0.5164228;
    a3=-2.737418;
    a4=16.07497;
    a5=-38.87607;
    a6=39.49944;
    a7=-14.17718;


    double A;
    A=D/(2*PI*PI);
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

double Lpeff_bouchiat(double pipette, double position, double &extension, double &z, double &force, double &forceTrap, double &Lpeff)
{
    double wlc; //no KTP coefficient
    double correction;

    z=(position-pipette)/L;
    extension=position-pipette;
    wlc=1.0/(4.0*pow((1.0-z),2))+z-0.25;
    correction=(a2*pow(z,2))+(a3*pow(z,3))+(a4*pow(z,4))+(a5*pow(z,5))+(a6*pow(z,6))+(a7*pow(z,7));
    force=KTP*(wlc+correction);
    Lpeff=(KT/force)*(wlc+correction);
    forceTrap=position*k;
}

void Langevin_position_WLC(double position, double force)
{
    position=position+(force*dts/G)+(D1*dts*position)+(sqrt(D2*dts)*rand_normal(0.0,1.0));
}

int main()
{
    double forceBuchiat,position,pipette,extension,forceTrap;
    double Lpeff;
    double z; //extension/L
    position=0.0;
    double time=0.0;

    initial_val();

    ofstream ofile;
    ofile.open("/Users/Naanaa/Desktop/buchiatnoise1.csv");
    do
    {
        for(int i=0;i<10;i++)
        {
            time=time+i*dts;
            Lpeff_bouchiat(pipette, position, extension, z, forceBuchiat, forceTrap, Lpeff);
            ofile<<pipette<<","<<position<<","<<extension<<","<<z<<","<<forceBuchiat<<","<<forceTrap<<","<<Lpeff<<endl;
            Langevin_position_WLC(position,forceBuchiat);
        }
        pipette=pipette+(0.1*1e-9);

    }while(time<Tmsr);

    ofile.close();
    return 0;
}

