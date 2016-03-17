/*this code is originally developed by David Boal and changes are made
by Naghmeh Rezaei in 2012

3D random walk: Monte Carlo method: Phantom chain (can intersect itself)
-----------------------------------------------------------------------*/


#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#define Lp 4.0 //Persistence length:flexural rigidity/kT=YI/kT

//constants used in random number function
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

#define PI 3.14159265

long int seed;

int seg; //number of joints in each chain
double bond_length; //separation between two joints (segment length);
int walks; //number of polymers (random walk chains to be generated)
double x[seg+1],y[seg+1],z[seg+1]; //position of each joint in the chain
int count[seg+1]; //counter for averaging in correlation calculations
double correlation[seg+1];
double Orient_diff[seg+1];

void initialize(void);
void random_walk(void);
void correlation_function(void);
float ran0(long *idum);


//initialize values
void initialize(void){

    printf ("Enter number of segments in eac chain: ");
    scanf ("%d" , &seg);

    printf ("Enter bond length: ");
    scanf ("%f" , &bond_length);

    printf ("Enter number of random walks (number of chains): ");
    scanf ("%d" , &walks);

    printf ("Enter persistence length (nm): ");
    scanf ("%f" , &Lp);

}//---------end of initilizer-----------



//uniform random number generator
float ran0(long *idum)
{

/* Minimal random number generator of Park and Miller, via numerical recipes
   has some correlations between successive numbers, not so great for sequential lattice
   calls  */
   long k;
   float ans;
   *idum ^= MASK;
   k=(*idum)/IQ;
   *idum = IA*(*idum-k*IQ)-IR*k;
   if(*idum < 0) *idum += IM;
   ans = AM*(*idum);
   *idum ^= MASK;
   return ans;

}//----------end of random numebr generator-----------



//random walk function generates each chain (polymer)
void random_walk(void)
{
   double cos_theta,sin_theta,phi; //theta: polar angle , phi: azimuthal angle
   double xn,yn,zn;
   double nn; //interaction between two nearest neighbours: dot product.

   x[0]=y[0]=z[0]=0.0;

   phi= 2.0*PI*ran0(&seed);
   cos_theta= 1.0-2.0*ran0(&seed);
   sin_theta= sqrt(1.0-cos_theta*cos_theta);

   x[1]=sin_theta*cos(phi);
   y[1]=sin_theta*sin(phi);
   z[1]=cos_theta;

   for(int i=2; i<=seg; i++)
   {
       for(int j=0; j<1000; j++) //1000 here play the role of a very big number
       {
           phi= 2.0*PI*ran0(&seed);
           cos_theta= 1.0-2.0*ran0(&seed);
           sin_theta= sqrt(1.0-cos_theta*cos_theta);

           xn=x[i-1]+sin_theta*cos(phi);
           yn=y[i-1]+sin_theta*sin(phi);
           zn=z[i-1]+cos_theta;

           nn=(xn-x[i-1])*(x[i-1]-x[i-2])
              +(yn-y[i-1])*(y[i-1]-y[i-2])
              +(zn-z[i-1])*(z[i-1]-z[i-2]);

           /* description for below if:
                (1-n.n)=1-cos(t=theta)
                =1-(1-t^2/2! + t^4/4! + ...)
                ~t^2/2
                so, Boltzmann probability to find system at
                bending energy Eb=kappa/2(t^2) is proportional to:
                exp(-Eb/Kt)*/

           if(exp(-Lp*(1.0-nn)) < ran0(&seed)) //metropolis algorithm
           {
               continue;
           }

           x[i]=xn;
           y[i]=yn;
           z[i]=zn;

           break;
          }
   }

} //----------end of walk -----------



//correlation function calculator
void correlation_function(void)
{
    int k;
    double nn; //t(x).t(x')
    double dn, dn2; //t(x)-t(x') , dn2=dn^2

    for (int i=1; i<seg/2; i++)
    {
        for (int j=0; j<seg/2; j++)
        {
            k=i+j;

            nn=(x[i]-x[i-1])*(x[k]-x[k-1])
                  + (y[i]-y[i-1])*(y[k]-y[k-1])
                  + (z[i]-z[i-1])*(z[k]-z[k-1]);

            count[j+1] +=1;
            correlation[j+1] +=nn;

            dn= (x[i]-x[i-1])-(x[k]-x[k-1]);
            dn2 = dn*dn;
            dn= (y[i]-y[i-1])-(y[k]-y[k-1]);
            dn2 += dn*dn;
            dn= (z[i]-z[i-1])-(z[k]-z[k-1]);
            dn2 += dn*dn;

            Orient_diff [j+1] +=dn2;
        }
    }

} //----------end of correlation function-----------



void main()
{
    double ree,ree2; //end-to-end distance of the polymer chain
    double avg_ree, avg_ree2;
    double var_ree;

    avg_ree=avg_ree2=0.0;

    for (int i=0; i<10; i++)
    {
        correlation[i]=0;
        count[i]=0;
        nnsq[i]=0;
    }

    printf ("\n\n\Analysis is performed on-line\n\n\n");

    //start loop over configuration sample
    for (int i=1; i<=walks; i++)
    {
        random_walk();
        correlation_function();

        ree2=(x[seg]*x[seg]+y[seg]*y[seg]+z[seg]*z[seg]);
        ree=sqrt(ree2);

        avg_ree2=avg_ree2+ree2;
        avg_ree=avg_ree+ree;

        if (walks%10==0)
        {
            printf ("%d random walks completed\n\n" , &i);
        }
    }

    printf ("%d random walks with %d segments\n\n" , &walks,seg);

    avg_ree=avg_ree/walks;
    avg_ree2=avg_ree2/walks;
    var_ree=avg_ree2-avg_ree*avg_ree;

    printf ("<ree>= %f \n", &avg_ree);
    printf ("<ree^2>= %f \n", &avg_ree2);
    printf (" variance of ree: <ree^2>-<ree>^2= %f \n", &var_ree);

    //calculating persistence length (Lp) from different results:

    printf ("\n Lp=kappa/kT= %f\n ", Lp);

    Lp= 0.5*(seg-sqrt(seg*seg-2.0*avg_ree2));
    printf (" <ree^2>~2LpLc-2Lp^2 ,so, Lp=0.5(Lc-sqrt(Lc^2-2<ree^2>))= %f\n" , Lp );

    Lp= 2.0*2.0
    printf (" <(t(s)-t(0))^2>=2s/Lp ,so, Lp=2s/<>= %f \n", Lp);

   printf(" del s   # in bin    <a(s)*a(0)>   <(a(s)-a(0))^2>\n");
   for(i=1; i<10; i++) {
      if(nubin[i] > 0) {cor[i]=cor[i]/nubin[i]; nnsq[i]=nnsq[i]/nubin[i];}
      printf("   %d       %d        %f       %f\n",i-1,nubin[i],cor[i],nnsq[i]);
   }

}
//-----------------------end of the code-------------------







