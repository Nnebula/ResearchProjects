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
double Tmsr,fs,dts,fnyq; //dt time step between computations
double k; //trap stiffness
double D,G,fc;
double D1,D2;

//functions:
void initial_val();
double rand_normal(double, double ); //generates random number with normal (gaussian) distribution
double find_extension(double, double, double);
void DNA_elastin_pulling_noise();


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

    //calculate constants:
    G=6*pi*1e-3*r; //all units SI
    fc=k/(2.0*pi*G);
    D=kb*temp/G;
    N=Tmsr*fs; //number of data generated

    dts=1.0/fs;
    fnyq=fs/2.0;

    D1=(-1.0)*2.0*pi*fc;//position is included in the main formula
    D2=2.0*D;

    double A;
    double tau;
    tau=(4*pi*pow(r,3)/(3*G));
    A=D/(2*pi*pi);
    cout<< "G:"<<G<<" , fc:"<<fc<<" , D:"<<D<<" , A:"<<A<<endl;
    cout<<"fs:"<<fs<<" , dts:"<<dts<<" , N:"<<N<<endl;
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

double find_extension(double myF,double contour_length, double persistence_length)
{
    double y_end,y_start;
    double y0;
    double guess_force;
    y_start=0.0;
    y_end=0.95*contour_length;

    y0=(y_end-y_start)/2;

    do{
        guess_force=(kb*temp/persistence_length)*((0.25/pow((1-(y0/contour_length)),2))+(y0/contour_length)-0.25);
        if(fabs(guess_force-myF)<0.001)
            return y0;

        if(myF>guess_force)
        {
            y_start=y0;
            y0=(y_end-y_start)/2;
        }

        if(myF<guess_force)
        {
            y_end=y0;
            y0=(y_end-y_start)/2;
        }

    } while(true);

}

//Trapped bead position with WLC force (pipette is moving with constant speed increasing the extension)
void DNA_elastin_pulling_noise()
{

    double L_DNA=340e-9;
    double P_DNA=50e-9;

    ofstream outWLC;
    outWLC.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/DNAelastinPullingNoise_20111215/1DNA340,50nm_100kHz_30pN_Vp50.txt");

    double time=0.0;

    double Xt; //trapped bead position
    Xt=0.0;
    double extension; //end to end distance of the DNA
    double WLC; // WLC force
    double force;
    int n=0;
    double Xp=0.0; // pipette bead position
//    double Xt_next;
    extension=0.0;

    do
    {
        force=k*Xt;
        extension=Xp-Xt;
        WLC=(kb*temp/P_DNA)*(1/(4*(pow(1-(extension/L_DNA),2)))+(extension/L_DNA)-0.25);

        outWLC<<n<<","<<time<<","<<Xt<<","<<Xp<<","<<extension<<","<<force<<","<<WLC<<endl;

        //next step:
        time=time+dts;
        Xt=Xt+D1*dts*Xt+(WLC*dts)/G+(sqrt(D2*dts)*rand_normal(0.0, 1.0));
        Xp=50e-9*time;
        n+=1;
    }while(force<30e-12);




//    WLC=(kb*temp/P_DNA)*(1/(4*(pow(1-(extension/L_DNA),2)))+(extension/L_DNA)-0.25);
//    Xt=WLC/k;
//    cout<< extension<<endl;

//    do
//    {
//        force=k*Xt;
//        WLC=force;

//        extension=find_extension(WLC,L_DNA,P_DNA);

//        Xt_next=Xt+D1*dts*Xt+(WLC*dts)/G+(sqrt(D2*dts)*rand_normal(0.0, 1.0));

//        outWLC<<n<<","<<time<<","<<Xt<<","<<Xp<<","<<extension<<","<<force<<","<<WLC<<endl;

//        Xp=100e-9*dts+Xp;
//        extension=Xp-Xt_next;
//        WLC=(kb*temp/P_DNA)*(1/(4*(pow(1-(extension/L_DNA),2)))+(extension/L_DNA)-0.25);
//        Xt=Xt_next;

//        n+=1;
//        time=time+dts;
//    }while(force<5e-12);

     outWLC.close();
}


///////////////////////////////////
int main()
{
    srand (time(NULL));

    initial_val(); //get initial values and calculates constants

    DNA_elastin_pulling_noise();

    return 0;

}

