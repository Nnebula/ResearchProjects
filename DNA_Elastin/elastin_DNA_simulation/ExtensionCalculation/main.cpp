//#include <QtCore/QCoreApplication>

#include <iostream>
#include <fstream>
#include <cmath>


const double kb=1.38e-23;
const int G=50000;

double P_DNA, P_elastin, L_DNA, L_elastin, T;

void Extension_calculations();

using namespace std;

int main()
{
    Extension_calculations();

    return 0;
}


void Extension_calculations()
{
    string fileName_DNA;
    string fileName_elastin;
    ofstream WLC;

    int DNA_number,elastin_number;

    cout << "enter the Temprature:"<< "  ";
        cin>> T;
        cout<< endl;

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

        cout << "enter file name(number of DNA/format:DNA:  -): ";
        cin >> fileName_DNA;


        cout << "enter file name(number of elastin/format:elastin:  ): ";
        cin >> fileName_elastin;

//        time_t rawtime;
//          struct tm * timeinfo;
//          char buffer [80];

//          time ( &rawtime );
//          timeinfo = localtime ( &rawtime );

//          strftime(buffer,80,"Now it's %I:%M%p.",timeinfo);
//        exit(0);

        WLC.open(("/Users/Naanaa/Documents/research/Trap_WLC_simulation/WLC_hybrid/"+ fileName_DNA+fileName_elastin +  ".txt").c_str());

        WLC << "DNA_extension"<<","<< "Force_DNA"<<","<<"elastin_extension"<<","<<"Force_DNA"<< ","<<"hybrid extension(nm): 1DNA" <<","<<"Force_hybrid_1DNA"<<","<<"hybrid extension(nm): 2DNA"<<","<<"Force_hybrid_2DNA"<<"\n";

    double myF[G],x_DNA[G], F_DNA, x_elastin[G], F_elastin;
    double dF;
    myF[0]=0;
    x_DNA[0]=x_elastin[0]=0;

    dF=0.1;
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
                if(fabs(myF[j]-F_DNA)<0.001)
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
                if(fabs(myF[i]-F_elastin)<0.0005)
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
            WLC << x_DNA[i]<<","<<myF[i]<<","<< x_elastin[i]<< ","<<myF[i]<<","<<x_DNA[i]+x_elastin[i]<<","<<myF[i]<<","<<(x_DNA[i]*DNA_number)+(x_elastin[i]*elastin_number)<<","<< myF[i]<< "\n";

            y1_elastin=0;
            y2_elastin=0.95*L_elastin;
            if(myF[i]>10)
            {
                break;
            }

        }

        WLC.close();


}
