#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

double Ncy[300000];
double y_e[300000];
double F[300000];

using namespace std;

int main() {
   ifstream fin;
   fin.open("/Users/Naanaa/Desktop/stdev_Analysis/Collagen_2-3-11_D3/2011-02-03-16-05-50_DAQ.txt", ios::in);
   string line;

   ofstream oFile;
   oFile.open("/Users/Naanaa/Desktop/stdev_Analysis/Collagen_2-3-11_D3/1_Stdev analysis_winSize20.txt");

   // skip line
   getline(fin, line);
   cout << "pooya2" << line << endl;

   // rest of the data
   int linenum=0;
   //int i;
   while (getline(fin, line)) {
      //getline(fin, line);
      //i++;

      //// replace comma with space
      //line.replace(line.find(","), 1, " ");
      istringstream linestream(line);

      // read tokens
      double a[6];

      linestream >> a[0] >> a[1] >> a[2] >> a[3] >> a[4] >> a[5];

      Ncy[linenum] = -5.0 * a[2] / a[3];
      y_e[linenum] = -5.08*a[5];
      y_e[linenum] = y_e[linenum]-Ncy[linenum];
      F[linenum] = Ncy[linenum]*112.93;
      linenum++;

//      Xe[i]=a;
//      Fe[i]=b;

//      oFile<< a << "," << b<< endl;

      //// skip everything small than 100.0
      //if (a<100.0) continue;

      // print
//      cout << "Row #" << linenum << ": " << a << "    " << b << endl;
   }
   cout << "Read " << linenum << " lines" << endl;



   const int winsize = 20;
   oFile<< "extension_average"<< " "<< "extension_stdev"<< " "<< "F_average"<< " "<< "Force_stdev"<< endl;

   for (int i=0; i<(linenum-winsize); i+=winsize) {
       double avgF = 0;
       double avgX = 0;
       for (int j=0; j<winsize; j++) {
           avgF += F[i+j];
           avgX += y_e[i+j];
       }
       avgF /= winsize;
       avgX /= winsize;

       double stdevF=0;
       double stdevX=0;
       for (int j=0; j<winsize; j++) {
           stdevF += pow(F[i+j]-avgF, 2);
           stdevX += pow(y_e[i+j]-avgX, 2);

       }
       stdevF = sqrt(stdevF/(winsize-1));
       stdevX = sqrt(stdevX/(winsize-1));

       oFile << avgX << "  "<<  stdevX+26.6322 << "  "<<avgF-7.94 << "  " << stdevF << endl;
   }


   oFile.close();

   return 0;

}
