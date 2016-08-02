#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

const int boxLen = 100;
double data[2][boxLen];

void getAvgDev(int col, int count, double &avg, double &stdev) {
   double sum = 0.0;
   for (int i=0; i<count; i++) {
      sum += data[col][i];
   }
   avg = sum / (float)count;

   sum = 0.0;
   for (int i=0; i<count; i++) {
      sum += pow(data[col][i] - avg, 2);
   }
   sum /= (float)count;
   stdev = sqrt(sum);
}

int main() {
   ifstream fin;
   fin.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/pullingDNA_20111027/pull_100speed_100kHz_elastin.txt", ios::in);
   string line;

   ofstream fout;
   fout.open("/Users/Naanaa/Documents/research/Trap_WLC_simulation/pullingDNA_20111027/SDEV_elastin.txt");

   int linenum = 0;
   while (getline(fin, line)) {

      // replace commas with space
      line.replace(line.find(","), 1, " ");
      int pos = line.find(",");
      while (pos!=string::npos) {
         line.replace(pos, 1, " ");
         pos = line.find(",", pos+1);
      }

      istringstream linestream(line);

      // read tokens
      double c1, c2, c3, c4, c5, c6, c7;
      linestream >> c1 >> c2 >> c3 >> c4 >> c5 >> c6 >> c7;

      // print
      // cout << "Row #" << linenum << ": " << c5 << " " << c6 << endl;

      // process
      data[0][linenum % boxLen] = c5;
      data[1][linenum % boxLen] = c6;
      if (linenum % boxLen==boxLen-1) {
         double d0avg, d0dev;
         double d1avg, d1dev;

         getAvgDev(0, boxLen, d0avg, d0dev);
         getAvgDev(1, boxLen, d1avg, d1dev);
         fout << d0avg << "," << d0dev << "," << d1avg << "," << d1dev << endl;
      }

      linenum++;
   }
   int remaining = linenum % boxLen;
   if (remaining>0) {
      double d0avg, d0dev;
      double d1avg, d1dev;

      getAvgDev(0, remaining, d0avg, d0dev);
      getAvgDev(1, remaining, d1avg, d1dev);
      fout << d0avg << "," << d0dev << "," << d1avg << "," << d1dev << endl;
   }
   fout.close();
}
