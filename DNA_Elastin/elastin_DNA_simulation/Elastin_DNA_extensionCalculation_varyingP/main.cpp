//#include <QtCore/QCoreApplication>

#include <iostream>
#include <fstream>
#include <cmath>


const double kb=1.38e-23;
const int G=50000;
const double T=298;

void DNA_calculate(double,double&); //for DNA calculation
void elastin_calculate(double , double , double , double&); // for peptide calculation

using namespace std;

int main()
{

    double P_elastin, L_elastin;
    ofstream hybrid;
    ofstream each;

    string fileName;
    cout<<"entre the file name:";
    cin>> fileName;

    cout<< "enter the contour lenth of peptide (nm):"<< "  ";
    cin>> L_elastin;
    cout<< endl;

    cout<< "enter the persistence lenth of peptide (nm):"<< "  ";
    cin>> P_elastin;
    cout<< endl;

    double x_DNA=0;
    double x_elastin=0;
    double myF=0;
    double df=0.0001;
    double dp=0.05;

    hybrid.open(("/Users/Naanaa/Desktop/ElastinDNA_22July2011/hybrid"+ fileName + ".txt").c_str());
    each.open(("/Users/Naanaa/Desktop/ElastinDNA_22July2011/single"+ fileName + ".txt").c_str());

    hybrid<<"Force";
    each<< "Force"<<","<<"x_1DNA";

    for(int j=0;j<20;j++)
    {
        hybrid<<","<<"Lp="<<P_elastin+dp*j;
        each<<","<<"Lp="<<P_elastin+dp*j;
    }

    hybrid<< endl;
    each<< endl;

    for(int j=0;j<G;j++)
    {
        myF=myF+df;

        DNA_calculate(myF, x_DNA);

        hybrid<<myF;
        each<< myF<<","<<x_DNA;

        for(double i=0;i<20;i++)
        {
          elastin_calculate(P_elastin+(i*dp), L_elastin, myF, x_elastin);
          hybrid<<","<<(2*x_DNA)+x_elastin;
          each<<","<<x_elastin;

        }


        hybrid<<endl;
        each<<endl;

        if(myF>10)
        {
            break;
        }
    }

    hybrid.close();
    each.close();

    return 0;
}


//DNA calculation

void DNA_calculate(double myF , double& y0_DNA)
{
    double y1_DNA=0;
    double y2_DNA=0.95*340;
    double F_DNA;

    for(int i=0;i<G;i++)
    {
       y0_DNA=(y1_DNA+y2_DNA)/2;

        do
        {
            F_DNA=(kb*T*1e21/50)*((1/(4*(1-(y0_DNA/340))*(1-(y0_DNA/340))))+(y0_DNA/340)-0.25);
            if(fabs(myF-F_DNA)<0.00005)
            {
                return;
            }

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

        }while (true);
    }
}


//elastin calculation

void elastin_calculate(double P, double L, double myF , double& y0_elastin)
{
    double y1_elastin=0;
    double y2_elastin=0.95*L;
    double F_elastin;

    for(int i=0;i<G;i++)
    {
       y0_elastin=(y1_elastin+y2_elastin)/2;

        do
        {
            F_elastin=(kb*T*1e21/P)*((1/(4*(1-(y0_elastin/L))*(1-(y0_elastin/L))))+(y0_elastin/L)-0.25);
            if(fabs(myF-F_elastin)<0.00005)
            {
                return;
            }

            if(myF>F_elastin)
            {
                y1_elastin=y0_elastin;
                y0_elastin=(y1_elastin+y2_elastin)/2;
            }

            if(myF<F_elastin)
            {
                y2_elastin=y0_elastin;
                y0_elastin=(y1_elastin+y2_elastin)/2;
            }

        }while (true);
    }
}






