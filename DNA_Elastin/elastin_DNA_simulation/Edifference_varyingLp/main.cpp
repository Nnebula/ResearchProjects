#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

//NOTE**) if the columns in the input file are seperated by space there is no need to go through lines 32:57
// in this case the file could be read directly by using: fin>>wave[i];

void file_read();


using namespace std;

int main() {

    file_read();

    return 0;
}

// function reads a file:
void file_read(){

   ifstream fin;
   fin.open("/Users/Naanaa/Desktop/ElastinDNA_22July2011/hybrid166,0.36,dp0.05.txt");

   int row_number=0;
   int col_num=0;

   ofstream oFile;
   oFile.open("/Users/Naanaa/Desktop/ElastinDNA_22July2011/Hybrid_Read_File_2.txt");

   //counting the number of columns:

   string line;
   getline(fin,line);
   int b=0;
   do{
       b=line.find(",",b+1);
       if(b== string::npos)
            break;
        col_num++;
    }while(true);
   cout<< col_num+1<<endl;


//   //put title for each column:
//   for(int i=1;i<col_num;i++)
//   {
//     oFile<<"L"<<0.36+(0.05*i)<<"-L"<<0.36<<","<<"L"<<0.36+(0.05*(i+1))<<"-L"<<0.36+(0.05*i)<<",";
//   }
//   oFile<<endl;


   // for all rows in the file:
   while(getline(fin,line))
   {
       row_number++;
       int search=0;
        // replace , with space
        do {
            search=line.find(",",search+1);
            if(search== string::npos)
                break;
            line.replace(search, 1, " ");

        } while(true);

        // read tokens

       double wave[col_num+1];
       double difference_L0[col_num+1];
       double difference_neighbor[col_num+1];
       istringstream read_line(line); //transforms line into an input string.

       for(int i=0;i<col_num;i++)
       {
           read_line>>wave[i];
       }


       for(int j=1;j<col_num;j++)
       {
           difference_L0[j]=wave[j+1]-wave[1]; //subtracts extensions of each column from column corresponding to L=0.36 (column#2) (Ex-E0.36)
           difference_neighbor[j]=wave[j+1]-wave[j];

           oFile<<difference_L0[j]<<","<<difference_neighbor[j]<<",";
       }
       oFile<<endl;

    }
   cout<< "row#: "<< row_number<< endl;


   oFile.close();

}
