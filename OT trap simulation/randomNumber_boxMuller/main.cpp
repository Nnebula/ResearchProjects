#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <iostream>
#include<fstream>

double rand_normal(double, double);

using namespace std;

int main()
{
    ofstream oFile;
    oFile.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/RandomNumber_20120122/1.txt");

    for(int i=0;i<1000000;i++)
    {
        oFile<<rand_normal(0.0,1.0)<<endl;

    }
   oFile.close();
   return 0;
}


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

