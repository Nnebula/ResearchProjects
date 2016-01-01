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

const double Tmsr=2;

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

    double L_DNA=680e-9;
    double P_DNA=50e-9;
    double L_elastin=166e-9;
    double P_elastin=0.36e-9;

    double extension_DNA,extension_elastin;
    double WLC_DNA, WLC_elastin;

    double extension_total; // extensionDNA+extensionElastin

    double x_trap; //position of the trapped bead
    double x_pipette;

    double time=0.0;
    double force; //force=x_trap*k

    double df; //force tollerance: f_elastin-f_DNA=df

    int m=0; //counts the number of searches

    double x_start,x_end,x_middle;

    ofstream outfile;
    outfile.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/BeadFluctuations@differentForces_20120203/DNA680_elastin166_extension0.00LDNA_100sec100k.txt");


    //initial values for x_pipette, x_trap, df
    df=0.001e-12;

    extension_DNA=0.00*L_DNA;

    WLC_DNA=(kb*temp/P_DNA)*(1/(4*(pow(1-(extension_DNA/L_DNA),2)))+(extension_DNA/L_DNA)-0.25);
    x_trap=WLC_DNA/k;

    cout<<"avg. force= "<<WLC_DNA<<endl;

    //find elastin extension at the above force:
    x_start=0.0;
    x_end=L_elastin;
    m=0;
    do{
        x_middle=(x_end-x_start)/2.0 + x_start;
        m+=1;

        extension_elastin=x_middle;

        WLC_elastin=(kb*temp/P_elastin)*(1/(4*(pow(1-(extension_elastin/L_elastin),2))) + (extension_elastin/L_elastin)-0.25);

        if(abs(WLC_elastin-WLC_DNA)<df){
            break;
        }

        if(WLC_DNA<WLC_elastin){
            x_end=x_middle;
        }

        if(WLC_DNA>WLC_elastin){
            x_start=x_middle;
        }

        if (m>5000) {
            cout << "FAILED!initially too many iterations for elastin calculations" << endl;
            exit(1);
        }

        if(extension_elastin>=L_elastin)
        {
            cout<<"WARNING in primary extension calculation!"<<endl<<"elastin EXTENSION OVER LIMIT!!"<<endl;
            exit(1);
        }

    }while(true);

    extension_total=extension_DNA+extension_elastin;
    x_pipette=x_trap+extension_total;

    cout<<"pipette position:   "<<x_pipette<<"          extension total:   "<<extension_total<<endl<<endl;



    outfile<<"time"<<","<<"x_trap"<<","<<"WLC_DNA"<<","<<"force"<<","<<"extension_total"<<endl;


    //calculating extension_DNA and extension_elastin for a given extension_total
    for(int n=0;n<10000000;n++)
    {

        if(n%1000000==0)
            cout<<"n= "<<n<<endl;

        x_start=0.0;
        x_end=extension_total;

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

        //write in the file
            outfile<<time<<","<<x_trap<<","<<WLC_DNA<<","<<force<<","<<extension_total<<endl;
            outfile.flush();

        //next time step
        time+=dts;

        x_trap=x_trap+D1*dts*x_trap+(WLC_DNA*dts)/G+(sqrt(D2*dts)*rand_normal(0.0, 1.0)); //langevin equation with thermal noise
        extension_total=x_pipette-x_trap;
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


