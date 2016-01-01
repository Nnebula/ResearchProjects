#include <iostream>
#include <fstream>
#include<cmath>

using namespace std;

int main()
{
    ofstream oFile;
    oFile.open("/Users/Naanaa/Desktop/elastin_DNA_simulation/plot_google_charts/chart.js");
    oFile << "function getData() { d = [ " << endl;

    for (int i=0; i<300; i++) {
        oFile << "[" << i << "," << 100.0*sin(i/50.0) << "]," << endl;
    }
    oFile << "]}" << endl;
    oFile.close();
}
