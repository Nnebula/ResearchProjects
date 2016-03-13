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

int main()
{
    double F,extension;
    extension=0.0;
    double KTP,F_DNA;
    double L=300*1e-9;
    KTP=kb*temp/(50*1e-9); //SI units

    ofstream ofile;
    ofile.open("/Users/Naanaa/Desktop/WLCpulling/WLCtest1.csv");

    for(int i=0;i<400;i++)
    {
        extension=extension+(i*(1e-9)/20);
        F=KTP*((1.0/(4.0*pow((1.0-(extension/L)),2)))+(extension/L)-0.25);
        F_DNA=(KTP)*((1/(4*(1-(extension/L))*(1-(extension/L))))+(extension/L)-0.25);
        ofile<<extension<<","<<F<<","<<F_DNA<<endl;
    }

    ofile.close();
    return 0;
}

