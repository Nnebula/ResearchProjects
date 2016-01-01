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
int avrg_length;

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

    avrg_length=166667;

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
    double WLC_DNA, WLC_elastin;
    double extension_elastin;

    double extension_total[avrg_length+1],extension_DNA[avrg_length+1]; // extensionDNA+extensionElastin

    double x_pipette,x_trap; //positions of the pipette and trapped beads
    double v_pipette; //pipette pulling speed

    double time=0.0;
    double force[avrg_length+1]; //force=x_trap*k

    double df; //force tollerance: f_elastin-f_DNA=df

    int n=0; //counter of total time steps
    int m=0; //counts the number of searches
    int l=1; //counts the number of data to be averaged over into a bin

    double sum_extension_total,sum_force,sum_extension_DNA;
    sum_extension_total=sum_force=sum_extension_DNA=0.0;
    double var_extension_total,var_force,var_extension_DNA,stdev_extension_total,stdev_force,stdev_extension_DNA;
    var_extension_total=var_force=var_extension_DNA=stdev_extension_DNA=0.0;
    double mean_extension_total,mean_force,mean_extension_DNA;

    double x_start,x_end,x_middle;

    ofstream outWLC;
    outWLC.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/DNAelastinPullingStdev_20120223/DNA680_elastin266_1usec_fs6Hz_v45.txt");

    outWLC<<"time"<<","<<"mean_extension_DNA"<<","<<"mean_extension_total"<<","<<"mean_force"<<","<<"stdev_extension_DNA"<<","<<"stdev_extension_total"<<","<<"stdev_force"<<
            ","<<"SEM_extension_DNA"<<","<<"SEM_extension_total"<<","<<"SEM_force"<<endl;

    //initial values for x_pipette, x_trap, v_pipette, df
    x_trap=0.0;
    x_pipette=60.0e-9;
    v_pipette=45.0e-9;
    df=0.001e-12;

    extension_total[1]=x_pipette-x_trap;

    //pulling starts here
    do{

        if(l<avrg_length)
        {
            n+=1;

            x_start=0.0;
            x_end=extension_total[l];

            //calculating extension_DNA and extension_elastin for a given extension_total
            m=0;
            while(true){
                x_middle=(x_end-x_start)/2.0 + x_start;
                m+=1;

                extension_DNA[l]=x_middle;
                WLC_DNA=(kb*temp/P_DNA)*(1/(4*(pow(1-(extension_DNA[l]/L_DNA),2)))+(extension_DNA[l]/L_DNA)-0.25);

                extension_elastin=extension_total[l]-extension_DNA[l];
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

                if (m>10000) {
                    cout << "FAILED! too many iterations" << endl;
                    exit(1);
                }
            }

            if(extension_elastin>=L_elastin || extension_DNA[l]>=L_DNA)
            {
                cout<<"WARNING!"<<endl<<"EXTENSION OVER LIMIT!!"<<endl;
                cout<< "elastin: "<< extension_elastin<< ","<< "DNA: "<<extension_DNA<<endl;
                exit(1);
            }

            //force calculation
            force[l]=k*x_trap;

//            //stiffness calculations:
//            k_DNA=(kb*temp/P_DNA)*(1/(2*L_DNA*(pow(1-(extension_DNA/L_DNA),3)))+1/L_DNA);
//            k_elastin=(kb*temp/P_elastin)*(1/(2*L_elastin*(pow(1-(extension_elastin/L_elastin),3)))+1/L_elastin);//derivative of WLC with respect to extension
//            k_hybrid=k_DNA*k_elastin/(k_DNA+k_elastin); // SI units


            sum_extension_total=sum_extension_total+extension_total[l];
            sum_force=sum_force+force[l];
            sum_extension_DNA=sum_extension_DNA+extension_DNA[l];

            l+=1;

            //next time step
            time+=dts;
            x_pipette=x_pipette+(v_pipette*dts);
            //x_trap=x_trap+D1*dts*x_trap+(WLC_DNA*dts)/G;  //no thermal noise included
            x_trap=x_trap+D1*dts*x_trap+(WLC_DNA*dts)/G+(sqrt(D2*dts)*rand_normal(0.0, 1.0)); //langevin equation with thermal noise
            extension_total[l]=x_pipette-x_trap;

            }

        //averaging and binnind each 10000 data points:
        else{

            mean_extension_total=sum_extension_total/avrg_length;
            mean_force=sum_force/avrg_length;
            mean_extension_DNA=sum_extension_DNA/avrg_length;

            for(int i=1;i<=avrg_length;i++)
            {
                var_extension_total=pow((extension_total[i]-mean_extension_total),2)+var_extension_total;
                var_force=pow((force[i]-mean_force),2)+var_force;
                var_extension_DNA=pow((extension_DNA[l]-mean_extension_DNA),2)+var_extension_DNA;
            }

            var_extension_total=var_extension_total/(avrg_length-1);
            var_force=var_force/(avrg_length-1);
            var_extension_DNA=var_extension_DNA/(avrg_length-1);

            stdev_extension_total=sqrt(var_extension_total);
            stdev_force=sqrt(var_force);
            stdev_extension_DNA=sqrt(var_extension_DNA);

            //writing in the file
            outWLC<<time<<","<<mean_extension_DNA*1e9<<","<<mean_extension_total*1e9<<","<<mean_force*1e12<<","<<stdev_extension_DNA*1e9<<","<<stdev_extension_total*1e9<<","<<stdev_force*1e12<<
                    ","<<stdev_extension_DNA*1e9/sqrt(avrg_length)<<","<<stdev_extension_total*1e9/sqrt(avrg_length)<<","<<stdev_force*1e12/sqrt(avrg_length)<<endl;
            outWLC.flush();

            //reset values
            l=1;
            sum_extension_total=sum_force=sum_extension_DNA=0.0;
            var_extension_total=var_force=var_extension_DNA=0.0;

        }


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

