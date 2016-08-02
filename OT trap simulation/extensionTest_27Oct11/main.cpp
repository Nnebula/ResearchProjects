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

int main()
{
    srand (time(NULL));
    int T=100;
    int set=1000;
    double position[set][T];
    double position2[set][T];
    double ave[T];

    ofstream myfile;
    myfile.open("/Users/Naanaa/Desktop/rand.csv");
    for(int j=0;j<set;j++)
    {
        for(int i=0;i<T;i++){
            position[j][i]=rand_normal(0.0,1.0) + (i>0 ? position[j][i-1] : 0.0);
            position2[j][i]=pow(position[j][i],2);
        }
    }
    for(int i=0;i<T;i++)
    {
        ave[i]=0.0;
        for(int j=0;j<set;j++)
        {
            ave[i]=position2[j][i]+ave[i];
        }
        ave[i]=ave[i]/set;
        myfile<<ave[i]<<endl;
    }
    myfile.close();
    return 0;
}
