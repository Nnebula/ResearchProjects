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
double fs,dts,fnyq; //dt time step between computations
double k,r; //trap stiffness
double D,G,fc;
double D1,D2;
double P_DNA,P_elastin,L_DNA,L_elastin;
double k_DNA,k_elastin,k_hybrid;

//functions:
void initial_val();
double rand_normal(double, double ); //generates random number with normal (gaussian) distribution
double find_extension(double, double, double);
void DNA_elastin_pulling();
void unit_conversion();


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


    cout<<"enter sampling frequency(Hz) fs:";
    cin>>fs;

    cout<<"enter trap stiffness(pN/um) k:";
    cin>>k;
    k=k*1e-6; //SI

    r=1e-6; //SI

    L_DNA=680e-9;
    P_DNA=50e-9;
    L_elastin=266e-9;
    P_elastin=0.36e-9;

    //calculate constants:
    G=6*pi*1e-3*r; //all units SI
    fc=k/(2.0*pi*G);
    D=kb*temp/G;

    dts=1.0/fs;
    fnyq=fs/2.0;

    D1=(-1.0)*2.0*pi*fc;//position is included in the main formula
    D2=2.0*D;

    double A;
    double tau;
    tau=(4*pi*pow(r,3)/(3*G));
    A=D/(2*pi*pi);
    cout<< "G:"<<G<<" , fc:"<<fc<<" , D:"<<D<<" , A:"<<A<<endl;
    cout<<"fs:"<<fs<<" , dts:"<<dts<<endl;
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
void DNA_elastin_pulling()
{

    double extension_DNA,extension_elastin;
    double WLC_DNA, WLC_elastin;

    double extension_total; // extensionDNA+extensionElastin

    double x_pipette,x_trap; //positions of the pipette and trapped beads
    double v_pipette; //pipette pulling speed

    double time=0.0;
    double force; //force=x_trap*k

    double f_old,WLC_old,extension_old;

    double df; //force tollerance: f_elastin-f_DNA=df

    int n=0; //counter of total time steps
    int m=0; //counts the number of searches

    double x_start,x_end,x_middle;

    ofstream outWLC;
    outWLC.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/DNAelastinPulling_20120222/DNA680_elastin266_1MHz_v45.txt");

    //initial values for x_pipette, x_trap, v_pipette, df
    x_trap=0.0;
    x_pipette=30.0e-9;
    v_pipette=45.0e-9;
    df=0.001e-12;

    force=x_trap*k;
    extension_total=x_pipette-x_trap;

    outWLC<<"time"<<","<<"extension_total"<<","<<"force"<<","<<"extension_DNA"<<","<<"extension_elastin"<<","<<"WLC_DNA"<<","<<"WLC_elastin"<<","<<"xDNA/L_DNA"<<","<<"xElastin/L_elastin"<<","<<"k_DNA"<<","<<"k_elastin"<<","<<"k_hybrid"<<endl;

    //pulling starts here
    do{
        f_old=force;
        WLC_old=WLC_DNA;
        extension_old=extension_DNA;

        n+=1;
        x_start=0.0;
        x_end=extension_total;

        //calculating extension_DNA and extension_elastin for a given extension_total
        m=0;
        while(true){
            x_middle=(x_end-x_start)/2.0 + x_start;
            m+=1;

            extension_DNA=x_middle;
            WLC_DNA=(kb*temp/P_DNA)*(1/(4*(pow(1-(extension_DNA/L_DNA),2)))+(extension_DNA/L_DNA)-0.25);

            extension_elastin=extension_total-extension_DNA;
            WLC_elastin=(kb*temp/P_elastin)*(1/(4*(pow(1-(extension_elastin/L_elastin),2))) + (extension_elastin/L_elastin)-0.25);

            if(abs(WLC_elastin-WLC_DNA)<df){
                break;
            }

            if(WLC_DNA>WLC_elastin){
                x_end=x_middle;
            }

            if(WLC_DNA<WLC_elastin){
                x_start=x_middle;
            }

            if (m>1000) {
                cout << "FAILED! too many iterations" << endl;
                exit(1);
            }
        }

        if(extension_elastin>=L_elastin || extension_DNA>=L_DNA)
        {
            cout<<"WARNING!"<<endl<<"EXTENSION OVER LIMIT!!"<<endl;
            cout<< "elastin: "<< extension_elastin<< ","<< "DNA: "<<extension_DNA<<endl;
            exit(1);
        }

        //force calculation
        force=k*x_trap;

        //stiffness calculations:
        k_DNA=(kb*temp/P_DNA)*(1/(2*L_DNA*(pow(1-(extension_DNA/L_DNA),3)))+1/L_DNA);
        k_elastin=(kb*temp/P_elastin)*(1/(2*L_elastin*(pow(1-(extension_elastin/L_elastin),3)))+1/L_elastin);//derivative of WLC with respect to extension
        k_hybrid=k_DNA*k_elastin/(k_DNA+k_elastin); // SI units

        //write in the file
        if(n%1000==0){
            outWLC<<time<<","<<extension_total*1e9<<","<<force*1e12<<","<<extension_DNA*1e9<<","<<extension_elastin*1e9<<","<<WLC_DNA*1e12<<","<<WLC_elastin*1e12<<","<<extension_DNA/L_DNA<<","<<extension_elastin/L_elastin<<","<<k_DNA*1e3<<","<<k_elastin*1e3<<","<<k_hybrid*1e3<<endl;
            outWLC.flush();
        }

        //next time step
        time+=dts;
        x_pipette=x_pipette+(v_pipette*dts);
        //x_trap=x_trap+D1*dts*x_trap+(WLC_DNA*dts)/G;  //no thermalnoise equation
        x_trap=x_trap+D1*dts*x_trap+(WLC_DNA*dts)/G+(sqrt(D2*dts)*rand_normal(0.0, 1.0)); //langevin equation with thermal noise
        extension_total=x_pipette-x_trap;


    }while(WLC_DNA<10e-12);

     outWLC.close();
}


///////////////////////////////////
int main()
{
    srand (time(NULL));

    initial_val(); //get initial values and calculates constants

    DNA_elastin_pulling();

    return 0;

}

