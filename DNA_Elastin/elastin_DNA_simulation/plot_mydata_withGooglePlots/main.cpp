#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

int main() {

    ifstream fin;
    fin.open("/Users/Naanaa/Desktop/elastin_DNA_simulation/Data/ElastinDNA_22July2011/hybrid166,0.36,dp0.05.txt");

    ofstream oFile;
    oFile.open("/Users/Naanaa/Desktop/plot_mydata_withGooglePlots/chart.js");
    oFile << "function getData() { d = [ " << endl;

    //"[" << i << "," << 100.0*sin(i/50.0) << "]," << endl;

    string line;
    getline(fin,line);

   while(getline(fin,line))
    {
       oFile<<"["<<line<<"],"<<endl;
   }


    oFile << "]}" << endl;
    oFile.close();

    return 0;
}
