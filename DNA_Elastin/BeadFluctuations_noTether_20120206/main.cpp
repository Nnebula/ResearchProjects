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

const double Tmsr=100;

//public variables:
double start_time = static_cast<double>(clock() / CLOCKS_PER_SEC);
double dts,fnyq; //dts time step between computations
double k; //trap stiffness
double D,G,fc;
double D1,D2;

//functions:
void initial_val();
double rand_normal(double, double ); //generates random number with normal (gaussian) distribution
void thermal_fluctuation();


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


    cout<<"enter time steps of the simulation (sec):";
    cin>>dts;

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

    D1=(-1.0)*2.0*pi*fc;//position is included in the main formula
    D2=2.0*D;

    cout<< "G:"<<G<<" , fc:"<<fc<<" , D:"<<D<<endl;
    cout<<"simulation time step (dts): "<<dts<<endl;
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

//Trapped bead position with WLC force (pipette is moving with constant speed increasing the extension)
void thermal_fluctuation()
{
    double x_trap=0.0; //position of the trapped bead

    double time=0.0;
    double trap_force; //force=x_trap*k
    double thermal_force; //random thermal force (stochastic term)

    double etha;
    etha=rand_normal(0.0, 1.0);

    ofstream outfile;
    outfile.open("/Users/Naanaa/Documents/Research/Data/simulation results/powerspectrum_notether/20131118.txt");


    outfile<<"time"<<","<<"x_trap"<<","<<"force"<<","<<"thermal_force_etha"<<endl;


    //calculating extension_DNA and extension_elastin for a given extension_total
    for(int n=0;n<10000000;n++)
    {

        if(n%1000000==0)
            cout<<"n= "<<n<<endl;

        //force calculation
        trap_force=k*x_trap;

        //write in the file
        outfile<<time<<","<<x_trap<<","<<trap_force<<","<<etha<<endl;
        outfile.flush();

        //next time step
        time+=dts;
        etha=rand_normal(0.0, 1.0); //thermal force
        thermal_force=(sqrt(G*kb*temp))*etha;

        x_trap=x_trap+D1*dts*x_trap+(sqrt(D2*dts)*etha); //langevin equation with no tether (no WLC force)
    }

    outfile.close();
}


///////////////////////////////////
int main()
{
    srand (time(NULL));

    initial_val(); //get initial values and calculates constants

    thermal_fluctuation();

    return 0;

}



