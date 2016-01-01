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
void WLC();


//calculates time in seconds
string current_time() {
    stringstream ss;
    double time = static_cast<double>(clock() / (double)CLOCKS_PER_SEC);
    ss.precision(3);
    ss << "[" << time-start_time << "] " ;
    return ss.str();
}


void WLC()
{

    double L_DNA=340e-9;
    double P_DNA=50e-9;
    double L_elastin=166e-9;
    double P_elastin=0.36e-9;

    double WLC_DNA,WLC_elastin;

    double dx=1e-9; //extesnion steps

    double extension_DNA,extension_elastin;
    extension_DNA=extension_elastin=-5.0e-9;

    ofstream outWLC;
    outWLC.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/DNAelastinPulling _20120123/WLC_DNA340_elastin166.txt");

    outWLC<<"extension_DNA"<<","<<"WLC_DNA"<<","<<"extension_elastin"<<","<<"WLC_elastin"<<endl;

    do{

        WLC_DNA=(kb*temp/P_DNA)*(1/(4*(pow(1-(extension_DNA/L_DNA),2)))+(extension_DNA/L_DNA)-0.25);

        WLC_elastin=(kb*temp/P_elastin)*(1/(4*(pow(1-(extension_elastin/L_elastin),2))) + (extension_elastin/L_elastin)-0.25);

        if(extension_DNA<340e-9)
        {
            outWLC<<extension_DNA<<","<<WLC_DNA;
        }

        if(extension_elastin<166e-9)
        {
            outWLC<<","<<extension_elastin<<","<<WLC_elastin;
        }

        outWLC<<endl;

        extension_DNA=extension_DNA+dx;
        extension_elastin=extension_elastin+dx;


    }while(WLC_DNA<100e-12 || WLC_elastin<100e-12);

     outWLC.close();
}


///////////////////////////////////
int main()
{
    srand (time(NULL));

    WLC();

    return 0;

}

