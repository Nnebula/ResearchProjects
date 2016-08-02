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
void Langevin_position_dicrete_Fokker(double,double);
void dft(double, double, double);

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


//langevin_discretized
void Langevin_position_dicrete_Fokker(double *time,double *position)
{
    time[0]=0.0;
    position[0]=0.0;

    for(int n=0;n<N-1;n++)
    {
        position[n+1]=position[n]+D1*dts*position[n]+(sqrt(D2*dts)*rand_normal(0.0, 1.0));
        time[n+1]=time[n]+dts;
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


//generates theoretical power spectrum (lorentzian)
void theoretical_PS(int k,double &freq,double &Pk_Lorentzian) //Pk theory
{
        freq=k/Tmsr;
        Pk_Lorentzian=D/(2*pi*pi*((fc*fc)+(freq*freq))); //lorentzian
}


//Trapped bead position with WLC force (pipette is moving with constant speed increasing the extension)
void Langevin_position_Fokker_WLC()
{
    ofstream outWLC;
    outWLC.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/Trap_WLC_2011102503/10MHz_0.1pN.txt");

    ofstream stdev;
    stdev.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/Trap_WLC_2011102503/10MHz_mean_0.1pN.txt");

    double time=0.0;

    double position; //trapped bead position
    double extension; //end to end distance of the DNA
    double WLC; // WLC force
    double force;

    double myF,F_DNA;
    double y1_DNA, y2_DNA, y0_DNA;
    y2_DNA=340e-9;
    y1_DNA=0.0;
    myF=0.1e-12;

    double mean_ex,mean_f,mean_WLC,mean_xt;
    mean_xt=mean_ex=mean_f=mean_WLC=0.0;

    y0_DNA=(y1_DNA+y2_DNA)/2;

    do{
        F_DNA=(kb*temp/50e-9)*(1/(4*(pow(1-(y0_DNA/340e-9),2)))+(y0_DNA/340e-9)-0.25);
        if(fabs(myF-F_DNA)<1e-15)
            break;

        if(myF>F_DNA)
        {
            y1_DNA=y0_DNA;
            y0_DNA=(y1_DNA+y2_DNA)/2;
        }

        if(myF<F_DNA)
        {
            y2_DNA=y0_DNA;
            y0_DNA=(y1_DNA+y2_DNA)/2;
        }

    } while(true);

    extension=y0_DNA;
    cout<< extension<<","<<F_DNA<<","<<myF<<endl;
    position=F_DNA/k;

    double position_last;

    for(int n=0;n<N;n++)
        {
            WLC=(kb*temp/50e-9)*(1/(4*(pow(1-(extension/340e-9),2)))+(extension/340e-9)-0.25);

            force=-position*k;

            outWLC<<n<<","<<time<<","<<position<<","<<extension<<","<<force<<","<<WLC<<endl;
            position_last=position;
            position=position+D1*dts*position+WLC*dts/G+(sqrt(D2*dts)*rand_normal(0.0, 1.0));
            extension=extension-(position-position_last);
            time=time+dts;
            mean_ex=mean_ex+extension;
            mean_f=mean_f+force;
            mean_WLC=mean_WLC+WLC;
            mean_xt=mean_xt+position;
        }
    stdev<<mean_xt/N<<","<<mean_ex/N<<","<<mean_f/N<<","<<mean_WLC/N;

    stdev.close();

    outWLC.close();

//    double stdev_ex,stdev_f,stdev_WLC;
//    double var_ex,var_f,var_WLC;
//    stdev_ex=stdev_f=stdev_WLC=0.0;
//    var_ex=var_f=var_WLC=0.0;
//    mean_ex=mean_f=mean_WLC=0.0;

//    for(int i=0;i<N-1;i++)
//    {
//        mean_ex=mean_ex+extension[i];
//        mean_f=mean_f+force[i];
//        mean_WLC=mean_WLC+WLC[i];
//    }
//    mean_ex=mean_ex/(N-1);
//    mean_f=mean_f/(N-1);
//    mean_WLC=mean_WLC/(N-1);

//    for(int i=0;i<N-1;i++)
//    {
//        var_ex=pow(mean_ex-extension[i],2)+var_ex;
//        var_f=pow(mean_f-force[i],2)+var_f;
//        var_WLC=pow(mean_WLC-WLC[i],2)+var_WLC;
//    }
//    var_WLC=var_WLC/(N-1);
//    var_ex=var_ex/(N-1);
//    var_f=var_f/(N-1);

//    stdev_ex=sqrt(var_ex);
//    stdev_f=sqrt(var_f);
//    stdev_WLC=sqrt(var_WLC);

//    stdev<<"var_extension"<<","<<"var_force"<<","<<"var_WLC"<<endl;
//    stdev<<var_ex<<","<<var_f<<","<<var_WLC<<endl<<endl;

//    stdev<<"stdev_extension"<<","<<"stdev_force"<<","<<"stdev_WLC"<<endl;
//    stdev<<stdev_ex<<","<<stdev_f<<","<<stdev_WLC<<endl;

//     stdev.close();
}


///////////////////////////////////
int main()
{
    srand (time(NULL));

//    ofstream outPS;
//    outPS.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/Trap_WLC_20111019/No_WLC.csv");

    initial_val(); //get initial values and calculates constants

//    double position_noWLC[N];
//    double time_noWLC[N];

//    Langevin_position_dicrete_Fokker(time_noWLC,position_noWLC); //calculates position of the trapped bead when there is no WLC force

//    double frequency[N/2];
//    double PS_noWLC[N/2];

//    dft(position_noWLC,frequency,PS_noWLC); //calculates the power spectrum of the trapped bead when there is no WLC force

    Langevin_position_Fokker_WLC();

//    for(int i=0;i<N;i++)
//    {
//        if(i<N/2)
//        {
//            outPS<<time_noWLC[i]<<","<<position_noWLC[i]<<","<<frequency[i]<<","<<PS_noWLC[i]<<endl;
//        }
//        else
//        {
//            outPS<<time_noWLC[i]<<","<<position_noWLC[i]<<endl;
//        }
//    }

//    outPS.close();

    return 0;
    }


