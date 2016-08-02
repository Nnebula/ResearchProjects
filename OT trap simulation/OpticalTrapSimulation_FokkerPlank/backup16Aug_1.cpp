#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include<fstream>
#include <time.h>

using namespace std;

int N;
double T,dt;
double k;
double D1,D2;
double D,G,fc;

const double pi=3.141592654;

void initial_val();
void position_cal(double ,double );
//double rand_num(int);
double rand_normal(double, double );
void dft(double, double , double,double);
void WLC(double, double);
void Langevin_WLC(double,double);


void initial_val()
{
    cout <<"enter measurement time(sec):";
    cin>>T;

    cout<<"enter dt(sec):";
    cin>>dt;

    cout<<"enter trap stiffness(pN/um):";
    cin>>k;
    k=k*1e-6; //SI

    double r;
    cout<<"enter bead radius(um):";
    cin>>r;
    r=r*1e-6; //SI

    //calculate constants:
    G=6*pi*1e-3*r; //all units SI
    fc=k/(2.0*pi*G);
    D=1.38*(1e-23)*298/G;
    N=T/dt;

    D1=(-1.0)*2.0*pi*fc;//position is included in the main formula
    D2=2.0*D;
    cout<< "G:"<<G<<" ,"<<"fc:"<<fc<<" ,"<<"D:"<<D<<" ,"<<"N:"<<N<<endl;
}

//random number generator
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

//solving langevin equation(finding position)

void position_cal(double *time,double *position)
{
    for(int n=0;n<N-1;n++)
    {
        time[n+1]=time[n]+dt;
        position[n+1]=position[n]+D1*dt*position[n]+(sqrt(D2*dt)*rand_normal(0.0, 2.0));
    }
}

//Discrete Fourier Transform
void dft(double position_in[], double *out_imag, double *out_real, double *k, double *mag2)
{
  long int i, j;
  double arg;
  double cosarg,sinarg;

  ofstream outPos;
  outPos.open("/Users/Naanaa/Documents/Programming/Research/OT trap simulation/OpticalTrapSimulation_FokkerPlank/data/PowerSpectrum_Aug16.csv");


  for(i=0; i<N; i++)
  {
    out_real[i] = 0.0;
    out_imag[i] = 0.0;
    arg = 2.0 * pi * (double)i / (double)N;
    for(j=0; j<N; j+=1)
    {
      cosarg = cos(j * arg);
      sinarg = sin(j * arg);
      out_imag[i] += (position_in[j] * sinarg);
      out_real[i] += (position_in[j] * cosarg);
    }
    out_imag[i]=out_imag[i]*dt;
    out_real[i]=out_real[i]*dt;
    k[i]=i*(1/T);
    mag2[i]=((out_real[i]*out_real[i])+(out_imag[i]*out_imag[i]))/T;
    outPos<<k[i]<<","<<mag2[i]<<endl;
  }
  outPos.close();
}

//WLC: extension from force calculations:
void WLC(double *extension, double *myF)
{
    string fileName_DNA;
    string fileName_elastin;
    ofstream WLC;

int DNA_number,elastin_number;

    cout<< "enter the contour lenth of DNA (nm):"<< "  ";
    cin>> L_DNA;
    cout<< endl;

    cout<< "enter the persistence lenth of DNA (nm):"<< "  ";
    cin>> P_DNA;
    cout<< endl;

    cout<< "enter the contour lenth of elastin (nm):"<< "  ";
    cin>> L_elastin;
    cout<< endl;

    cout<< "enter the persistence lenth of elastin (nm):"<< "  ";
    cin>> P_elastin;
    cout<< endl;

    cout << "enter the number of DNA in hybrid system"<< "  ";
    cin>> DNA_number;
    cout<< endl;

    cout << "enter the number of elastinin hybrid system"<< "  ";
    cin>> elastin_number;
    cout<< endl;

    cout << "enter file name(number of DNA)format:  DNA:  -: ";
    cin >> fileName_DNA;

    cout << "enter file name(number of elastin/format:elastin:  ): ";
    cin >> fileName_elastin;

    WLC.open(("/Users/Naanaa/Desktop/ElastinDNA_19July2011/ContourLength_Constant/"+ fileName_DNA+fileName_elastin +  ".txt").c_str());

    WLC << "DNA_extension"<<","<< "Force_DNA"<<","<<"elastin_extension"<<","<<"Force_DNA"<< ","<<"hybrid extension(nm): 1DNA" <<","<<"Force_hybrid_1DNA"<<","<<"hybrid extension(nm): 2DNA"<<","<<"Force_hybrid_2DNA"<<"\n";

    double x_DNA[G], F_DNA, x_elastin[G], F_elastin;
    double dF;
    myF[0]=0;
    x_DNA[0]=x_elastin[0]=0;

    dF=0.00010;
    double y1_DNA, y1_elastin, y2_DNA, y2_elastin,y0_DNA, y0_elastin;
    y1_DNA=y1_elastin=0;
    y2_DNA=0.95*L_DNA;
    y2_elastin=0.95*L_elastin;

    //DNA calculations:
    for(int j=1; j<G; j++)
    {
            myF[j]=myF[j-1]+dF;
            y0_DNA=(y1_DNA+y2_DNA)/2;

            do{
                F_DNA=(kb*T*1e21/P_DNA)*((1/(4*(1-(y0_DNA/L_DNA))*(1-(y0_DNA/L_DNA))))+(y0_DNA/L_DNA)-0.25);
                if(fabs(myF[j]-F_DNA)<0.00005)
                    break;

                if(myF[j]>F_DNA)
                {
                    y1_DNA=y0_DNA;
                    y0_DNA=(y1_DNA+y2_DNA)/2;
                }

                if(myF[j]<F_DNA)
                {
                    y2_DNA=y0_DNA;
                    y0_DNA=(y1_DNA+y2_DNA)/2;
                }
            } while(true);

            x_DNA[j]=y0_DNA;

            y1_DNA=0;
            y2_DNA=0.95*L_DNA;
            if(myF[j]>10)
            {
                break;
            }
        }

    //elastin calculations:
    for(int i=1; i<G; i++)
    {
            y0_elastin=(y1_elastin+y2_elastin)/2;

            do{
                F_elastin=(kb*T*1e21/P_elastin)*((1/(4*(1-(y0_elastin/L_elastin))*(1-(y0_elastin/L_elastin))))+(y0_elastin/L_elastin)-0.25);
                if(fabs(myF[i]-F_elastin)<0.00005)
                    break;

                if(myF[i]>F_elastin)
                {
                    y1_elastin=y0_elastin;
                    y0_elastin=(y1_elastin+y2_elastin)/2;
                }

                if(myF[i]<F_elastin)
                {
                    y2_elastin=y0_elastin;
                    y0_elastin=(y1_elastin+y2_elastin)/2;
                }
            }while(true);

            x_elastin[i]=y0_elastin;
            extension[i]=x_elastin[i]+x_DNA[i];

  //          WLC << x_DNA[i]<<","<<myF[i]<<","<< x_elastin[i]<< ","<<myF[i]<<","<<x_DNA[i]+x_elastin[i]<<","<<myF[i]<<","<<(x_DNA[i]*DNA_number)+(x_elastin[i]*elastin_number)<<","<< myF[i]<< "\n";

            y1_elastin=0;
            y2_elastin=0.95*L_elastin;
            if(myF[i]>10)
            {
                break;
            }
        }

        WLC.close();
}

//solving langevin equation with WLC force(finding position)

void Langevin_WLC(double *time,double *position)
{
    for(int n=0;n<N-1;n++)
    {
        time[n+1]=time[n]+dt;
        position[n+1]=position[n]+D1*dt*position[n]+(sqrt(D2*dt)*rand_normal(0.0, 2.0));
    }
}

////////////
int main()
{
    srand (time(NULL));

    initial_val(); //get initial values and calculates constants

    double position[N];
    double time[N];
    position[0]=0;
    time[0]=0.0;
    double FT_pos_real[N],FT_pos_imag[N];
    double Pk[N];
    double f[N];

    position_cal(time,position);//Langevinequation solution (fokker_plank)

//    dft(position, FT_pos_real, FT_pos_imag, f, Pk); //Discrete fourier transform

    WLC(extension, myF);//WLC: extension from Force;

    return 0;
}
