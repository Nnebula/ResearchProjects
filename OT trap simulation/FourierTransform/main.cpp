#include <cmath>
#include <cstdlib>
#include<iostream>
#include<fstream>

using namespace std;

//#define TRUE 1
//#define FALSE 0

// Take, eat, this is my code which is hacked for you...
// (In case you didn't grasp that one consider this public domain)
// nolandda Thu Mar 11 03:49:18 EST 2004

/*
  Discrete Fourier Transform
*/
//out_imag & out_real: fourier transformed
void dft(long int length, double input_sample[], double *out_imag, double *out_real)
{
  long int i, j;
  double arg;
  double cosarg,sinarg;

  for(i=0; i<length; i++)
  {
    out_real[i] = 0;
    out_imag[i] = 0;
    arg = -1.0 * 2.0 * 3.141592654 * (double)i / (double)length;
    for(j=0; j<length; j+=1)
    {
      cosarg = cos(j * arg);
      sinarg = sin(j * arg);
      out_imag[i] += (input_sample[j] * sinarg);
      out_real[i] += (input_sample[j] * cosarg);
//      temp_imag[i] += (real_sample[j] * sinarg + imag_sample[j] * cosarg);
//      temp_real[i] += (real_sample[j] * cosarg - imag_sample[j] * sinarg);
    }
  }

}


///*
//  Inverse Discrete Fourier Transform
//*/

//int inverse_dft(long int length, double real_sample[], double imag_sample[])
//{
//  long int i, j;
//  double arg;
//  double cosarg,sinarg;
//  double *temp_real=NULL,*temp_imag=NULL;

//  temp_real = calloc(length, sizeof(double));
//  temp_imag = calloc(length, sizeof(double));
//  if (temp_real == NULL || temp_imag == NULL)
//  {
//    return(FALSE);
//  }

//  for(i=0; i<length; i+=1)
//  {
//    temp_real[i] = 0;
//    temp_imag[i] = 0;
//    arg = 2.0 * 3.141592654 * (double)i / (double)length;
//    for(j=0; j<length; j+=1)
//    {
//      cosarg = cos(j * arg);
//      sinarg = sin(j * arg);
//      temp_real[i] += (real_sample[j] * cosarg - imag_sample[j] * sinarg);
//      temp_imag[i] += (real_sample[j] * sinarg + imag_sample[j] * cosarg);
//    }
//  }

//  /* Copy the data back */
//  for (i=0; i<length; i+=1)
//  {
//    real_sample[i] = temp_real[i] / (double)length;
//    imag_sample[i] = temp_imag[i] / (double)length;
//  }

//  free(temp_real);
//  free(temp_imag);
//  return(TRUE);
//}

int main()
{
//  double* dr = calloc(100, sizeof(double));
//  double* di = calloc(100, sizeof(double));
//  dft(100, dr, di);
//  inverse_dft(100, dr, di);
    ofstream FT_file;
    FT_file.open("/Users/Naanaa/Desktop/OT_simulations/FT.csv");

    int length = 100;
    double out_real[length];
    double out_imag[length];
    double d[length];
    for(int i=0;i<length;i++)
    {
        d[i]=sin(2*3.14*i*5/length)*10;
    }
    dft(length,d,out_real,out_imag);

    for(int i=0;i<length;i++)
    {
        FT_file<<i<<","<<d[i]<<","<<out_real[i]<<","<<out_imag[i]<<endl;
    }
    FT_file.close();

  return 0;
}
