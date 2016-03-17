#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* 3D chain under stress */
/* MONTE CARLO */

#define BK 4.0  // bending constant of chain (kappa/kT = Lp)
#define TMAX 2.0  // TMAX is max tether length squared
// TMAX = 4 for 2D in general,  TMX = 2 for 3D
#define NSEG 20  // number of segments
// segment size is 1.1 times hard core diameter
#define STRESS 0.0  // force for chains, pressure for rings
#define NKEEP 300   // number of configs in analysis
#define NRELAX 210  // number of initial configs discarded
#define SWEEP 1000  // number of sweeps between output
#define IFANA 1    // set = 1 for on-the-fly analysis and no output
#define ipn (NSEG+1) //number of particles (joints)
//numbering of the particles start from 0 (first) to NSEG(last)

#define PI 3.14159265

// random generator function constants
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

long seed;

int ntot= NRELAX + NKEEP; //total number of configurations

double stres;

// variables used in main and analysis;
int config; // counter
double avg_rg2; // avg of radius of gyration for all KEEP chain configs
double avg_ree2; // avg of end to end distance over all KEEP chain configs
double rpos;

double x[ipn], y[ipn], z[ipn]; // joint positions

double ncoor[ipn]; // number of connections to each particle (1 at ends and 2 for middle)
double neighbor[ipn][2]; // particle number (i) of previous and next neighbor
//[i][0] previous neighbor & [i][1] next neighbor for ith particle
int nloc[ipn]; // numbering of particles occupying hardcore region of each particle
int nocc[ipn][50];

FILE *out_SA1;

void initialize(void);
void lookup(void);
void position(void);
void analysis(void);
void output(void);
float ran0(long *idum);


float ran0(long *idum) {
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
}// -----------------end of ran0--------------------


void initialize(void) { // setting up initial configuration of the chain as an open circle
    double step, radius;
    double ang; // angle of each particle on a circle (ring)
    int i;

    seed=234243;

    step=2.0*PI/ipn;
    radius=0.5*1.2*ipn/PI;
    ang= -0.5*PI + 0.5*step;

    for ( i=0; i<ipn; i++) { // setting up initial position of each particle on a ring
        ang += step;

        x[i]= radius* cos(ang);
        y[i]= radius* sin(ang);
        z[i]= 0.0;
    }

    for ( i=0; i<ipn; i++) {
        if (i==0 || i==(ipn-1) ) {
            ncoor[i]=1;
        }
        ncoor[i]=2;
        neighbor[i][0]= i-1; // previous neighbor of i
        neighbor[i][1]= i+1; // next neighbor of i
    }
}  // --------------------end of initialize----------------------


void lookup(void) {

    int i,j;
    double rmax2=4.0;
    double xcm, ycm, zcm;
    xcm= ycm= zcm=0.0;

    int noci, nocj; // numbering of particles occupying hardcore region of particle i, or j
    double dx,dy,dz; // distance between particle i and j
    double spacing2;

    for ( i=0; i<ipn; i++) { // zero the centre of mass positions
        nloc [i]=0.0;

        xcm +=x[i];
        ycm +=y[i];
        zcm +=z[i];
    }

    xcm /= ipn;
    ycm /= ipn;
    zcm /= ipn;

    for ( i=0; i<ipn; i++) {
        x[i] +=-xcm;
        y[i] +=-ycm;
        z[i] +=-zcm;
    }

    for ( i=0; i<ipn; i++) {
        for ( j=0; j<i; j++) {
            dx= x[i]-x[j];
            dy= y[i]-y[j];
            dz= z[i]-z[j];

            spacing2= dx*dx+ dy*dy+ dz*dz;

            if (spacing2 > rmax2) {
                continue;
            }
            else {
                noci= nloc[i];
                nocc[i][noci]=j;

                nocj= nloc[j];
                nocc[j][nocj]=i;

                nloc[i]+=1;
                nloc[j]+=1;
            }
        }
    }
}// --------------------end of lookup---------------------


void position(void) {

    int i, j, k;
    int iflag;
    int nn; //particle number which is close to particle i (e.g. previous=i-1 & next=i+1)
    double temp_pos[3]; // temporary(trial) position of each particle (0=x, 1=y, 2=z)
    double dx,dy,dz; //distance from previous and next neighbors
    double spacing2;
    double ds= 0.05; //displacement of each particle in each step from -ds/2 to +ds/2

    double dw; // work done by the last segment on the applied force (stress)

    // trial moves on particle positions
    for ( i=1; i<ipn; i++) {
        iflag= 0;

        temp_pos[0]= x[i]+ds*(0.5-(ran0(&seed))); //x
        temp_pos[1]= y[i]+ds*(0.5-(ran0(&seed))); //y
        temp_pos[2]= z[i]+ds*(0.5-(ran0(&seed))); //z

        // test for maximum tether length

        if (i==ipn-1) {
            dx= temp_pos[0]-x[ipn-2];
            dy= temp_pos[1]-y[ipn-2];
            dz= temp_pos[2]-z[ipn-2];

            spacing2= dx*dx+ dy*dy+ dz*dz;
            if (spacing2> TMAX) {
                iflag=1;
            }
        }
        else {
            for ( j=0; j<ncoor[i]; j++) {
                nn=neighbor[i][j];
                dx= temp_pos[0]-x[nn];
                dy= temp_pos[1]-y[nn];
                dz= temp_pos[2]-z[nn];

                spacing2= dx*dx+ dy*dy+ dz*dz;

                if (spacing2> TMAX) {
                    iflag=1;
                    break;
                }
            }  //end of j
        }
        if (iflag ==1) {
            continue;
        }
        // test for hardcore region
        for ( k=0; k<nloc[i]; k++) {
            nn= nocc[i][k];
            dx= x[nn]-temp_pos[0];
            dy= y[nn]-temp_pos[1];
            dz= z[nn]-temp_pos[2];

            spacing2= dx*dx+ dy*dy+ dz*dz;

            if (spacing2<1.0) {
                iflag=1;
                break;
            }
        }
        if (iflag== 1) {
            continue;
        }
        //test force on the chain
        if (i== (ipn-1)) {
            dw= stres*(temp_pos[i]-x[i]); //work on the particle by the force (stress)

            if (dw < 0.0) {
                if (exp(-dw) > ran0(&seed)) {
                    iflag=1;
                }
            }
        }

        if (iflag ==1) {
            continue;
        }
        // acepting the move
        x[i]= temp_pos[0];
        y[i]= temp_pos[1];
        z[i]= temp_pos[2];
    }// end of i loop over particles
}//--------------end of position-----------------


int main()
{
    int j;
    out_SA1=fopen("/Users/Naanaa/Documents/Programming/ResearchCodes/MonteCarlo/DATA/SA_f0_20NSEG_ds005_sweep1000.csv","w");

    initialize();
    stres=STRESS;

    avg_ree2=avg_rg2=0.0;

    // start loop over configuration sample
    for (config=0; config<ntot; config++) {
        rpos= 0.0;
        for (j=0; j<SWEEP; j++) {
            if (j%10==0) {
                lookup();
            }
            position();
        } // end of j
        if (config < NRELAX) {
            continue;
        }
        analysis();
    } // end of "config" loop over configurations

    fclose(out_SA1);
}// ---------------end of main routine----------------


void analysis(void) {

    int i;
    double ree2;
    double rg2;
    double xcm, ycm, zcm;

    // zero the centre of mass of the chain
    xcm= ycm= zcm= 0.0;
    for ( i=0; i<ipn; i++) {
        xcm += x[i];
        ycm += y[i];
        zcm += z[i];
    }
    xcm /= ipn;
    ycm /= ipn;
    zcm /= ipn;

    // sum of (end to end distance^2) over KEEP configs
    ree2 = (x[NSEG]-x[0])*(x[NSEG]-x[0])
           + (y[NSEG]-y[0])*(y[NSEG]-y[0])
           + (z[NSEG]-z[0])*(z[NSEG]-z[0]);
    avg_ree2 +=ree2;

    // average of (radius of gyration^2) for each chain
    for ( i=0; i<ipn; i++) {
        rg2= (x[i]-xcm)*(x[i]-xcm)
             + (y[i]-ycm)*(y[i]-ycm)
             + (z[i]-zcm)*(z[i]-zcm);
    }
    rg2 /= ipn;
    avg_rg2 += rg2;

    for ( i=0; i<ipn; i++) {
        fprintf (out_SA1, "%9.6f , %9.6f , %9.6f , ", x[i], y[i], z[i]);
    }
    fprintf (out_SA1, "\n");
    fflush (out_SA1);

    if (config == (ntot-1) ) {
        avg_ree2 /= NKEEP;
        avg_rg2 /= NKEEP;
    }
    if(config == (ntot-1) ) {
        avg_ree2 /= NKEEP;
        avg_rg2 /= NKEEP;
        printf("   <ree2>=%9.7f   <rg2>=%9.7f\n", avg_ree2, avg_rg2);
    }
} //-----------------end of analysis--------------------


//void output(void)
//{
//   int i, nsegm1;
//   double ree2;
///* some simple analysis */
//   nsegm1=NSEG-1;
//   ree2 = (x[nsegm1]-x[0])*(x[nsegm1]-x[0]) + (y[nsegm1]-y[0])*(y[nsegm1]-y[0])
//      + (z[nsegm1]-z[0])*(z[nsegm1]-z[0]);
///* write out the configurations */
//   fprintf(fpo," stres = %5.1f\n",stres);
//   fprintf(fpo,"%d   %d\n",NSEG,ipn);
//   for(i=0; i<ipn; i++) fprintf(fpo,"%f %f %f\n", x[i],y[i],z[i]);
//   fflush(fpo);
//   rpos=rpos/(SWEEP*ipn);
//   fprintf(fpl,
//     " # %d  rpos=%6.3f  stres=%5.1f   ree2=%6.2f\n",config,rpos,stres,ree2);
//   fflush(fpl);
//}

