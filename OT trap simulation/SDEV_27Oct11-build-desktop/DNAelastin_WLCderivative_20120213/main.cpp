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

//public variables:
double kbT;
double kb;
double temp;

double L_DNA;
double P_DNA;
double L_elastin;
double P_elastin;

double force_max,df,force_step;



//functions:
void DNA_elastin_stiffness();
void theory(); //stiffness at very small forces predicted by formula: 3kbT/2PL
void unit_conversion(); //converts units from SI(current unit to pN, nm, um )

void initial_values()
{
    kb=1.38*1e-23;//SI unit
    temp=298; //(kelvin)

    L_DNA=680e-9;
    P_DNA=50e-9;
    L_elastin=166e-9;
    P_elastin=0.36e-9;

    kbT=kb*temp;

    force_max=12e-12;
    df=0.00001e-12;
    force_step=0.001e-12;

}


void unit_conversion(){
    kbT=kb*temp*1e21;

    L_DNA=L_DNA*1e9;
    L_elastin=L_elastin*1e9;
    P_DNA=P_DNA*1e9;
    P_elastin=P_elastin*1e9;

    force_step=force_step*1e12;
    df=df*1e12;
    force_max=force_max*1e12;

}


void theory()
{
    double k_DNA,k_elastin,k_hybrid;

    k_DNA=(3*kbT)/(2*P_DNA*L_DNA);
    k_elastin=(3*kbT)/(2*P_elastin*L_elastin);
    k_hybrid=(k_DNA*k_elastin)/(k_elastin+k_DNA);

    cout<< "k_elastin:  "<<k_elastin<<endl;
    cout<< "k_DNA:  "<<k_DNA<<endl;
    cout<< "k_hybrid:  "<<k_hybrid<<endl;

    cout<<"kbT:  "<<kbT<<endl;
}

////////////////////////////
void DNA_elastin_stiffness()
{

    double extension_DNA,extension_elastin;
    double WLC_DNA, WLC_elastin;
    double k_elastin,k_DNA,k_hybrid;
    double force=0.0;

    int m=0; //counts the number of searches

    double x_start,x_end,x_middle;

    ofstream outWLC;
    outWLC.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/DNAelastin_stiffness_20120213/leastin166_DNA680_SI.txt");

    force;
    outWLC<<"extension_elastin"<<","<<"extension_DNA"<<","<<"extension_total"<<","<<"force"<<","<<"k_elastin"<<","<<"k_DNA"<<","<<"k_hybrid"<<endl;

    do{
        force+=force_step;

        //calculating extension_DNA at given forces
        if(extension_DNA<L_DNA)
        {
            m=0;
            x_start=0.0;
            x_end=L_DNA;

            while(true){
                x_middle=(x_end-x_start)/2.0 + x_start;
                m+=1;

                extension_DNA=x_middle;
                WLC_DNA=(kbT/P_DNA)*(1/(4*(pow(1-(extension_DNA/L_DNA),2)))+(extension_DNA/L_DNA)-0.25);

                if(abs(force-WLC_DNA)<df){
                    break;
                }

                if(WLC_DNA>force){
                    x_end=x_middle;
                }

                if(WLC_DNA<force){
                    x_start=x_middle;
                }

                if (m>1000) {
                    cout << "DNA FAILED! too many iterations" << endl;
                    exit(1);
                }
            }
            k_DNA=(kbT/P_DNA)*(1/(2*L_DNA*(pow(1-(extension_DNA/L_DNA),3)))+1/L_DNA); //derivative of WLC with respect to extension

        }

        //calculating extension_elastin at given forces
        if(extension_elastin<L_elastin)
        {
            m=0;
            x_start=0.0;
            x_end=L_elastin;

            while(true){
                x_middle=(x_end-x_start)/2.0 + x_start;
                m+=1;

                extension_elastin=x_middle;
                WLC_elastin=(kbT/P_elastin)*(1/(4*(pow(1-(extension_elastin/L_elastin),2))) + (extension_elastin/L_elastin)-0.25);

                if(abs(force-WLC_elastin)<df){
                    break;
                }

                if(WLC_elastin>force){
                    x_end=x_middle;
                }

                if(WLC_elastin<force){
                    x_start=x_middle;
                }

                if (m>1000) {
                    cout << "elastin FAILED! too many iterations" << endl;
                    exit(1);
                }
            }
            k_elastin=(kbT/P_elastin)*(1/(2*L_elastin*(pow(1-(extension_elastin/L_elastin),3)))+1/L_elastin); //derivative of WLC with respect to extension
        }

        k_hybrid=(k_DNA*k_elastin)/(k_DNA+k_elastin);

        outWLC<<extension_elastin<<","<<extension_DNA<<","<<extension_DNA+extension_elastin<<","<<force<<","<<k_elastin<<","<<k_DNA<<","<<k_hybrid<<endl;


    }while(force<force_max);

     outWLC.close();
}


///////////////////////////////////
int main()
{
    initial_values();
   // unit_conversion();

    theory();

    DNA_elastin_stiffness();

    return 0;

}


