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

int main()
{
    double force,extension;
    extension=0.0;
    double L=3.958*1e-6; //SI units
    double KTP=KB*TEMP/(50*1e-9); //SI units
    double KT=KB*TEMP; //SI units
    double Lpeff;
    double wlc; //no KTP coefficient
    double a2,a3,a4,a5,a6,a7;

    a2=-0.5164228;
    a3=-2.737418;
    a4=16.07497;
    a5=-38.87607;
    a6=39.49944;
    a7=-14.17718;

    cout<< a2+a3+a4+a5+a6+a7;

    double z; //extension/L


    ofstream ofile;
    ofile.open("/Users/Naanaa/Documents/Programming/Research/WLC_Lpeff/WLC_Lpeff_buchiat/Data/3.958umDNA.csv");

    ofile<< "extension(um)" <<","<< "extension/L"<<","<<  "force(pN)" <<","<< "Lpeff(nm)" <<endl;

    do
    {
        z=extension/L;
        wlc=1.0/(4.0*pow((1.0-z),2))+z-0.25;
        force=KTP*(wlc+(a2*pow(z,2))+(a3*pow(z,3))+(a4*pow(z,4))+(a5*pow(z,5))+(a6*pow(z,6))+(a7*pow(z,7)));
        Lpeff=(KT/force)*(wlc+(a2*pow(z,2))+(a3*pow(z,3))+(a4*pow(z,4))+(a5*pow(z,5))+(a6*pow(z,6))+(a7*pow(z,7)));

        ofile<< extension*1e6 <<","<< z<<","<<  force*1e12 <<","<< Lpeff*1e9 <<endl;
        extension=extension+(0.05*(1e-9));
    }while(force<=10*1e-12);

    ofile.close();
    return 0;
}

