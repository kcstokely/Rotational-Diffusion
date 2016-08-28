/******************************************************************************

  ROTATIONAL DIFFUSION:

    This program creates a file listing the {x, y, z} coordinates
    of a point on the surface of a sphere as a function of time.
    
    Tau values are in terms of the rank-2 correlation functions.

******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  C      2.35482005  // full width at half max by stnd dev
#define  PI     3.14159265  // circumference by double the radius
#define  DEC    16          // number of decimal places in output
#define  AXIS   0           // set to 0 for correct rotation axis

double          p[4];
extern void     rotate(double);
extern double   sample_rayleigh();
extern double   sample_gaussian();

extern int      check_flags(int, char**);

int             print_tau,print_head;
int             seed,runs,start,trajectory;
double          tau_median,tau_fwhm,tau_corr;
double          xch_median,xch_fwhm,xch_corr;

///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{

  /////////////////////////////////////////////////////////////////////////////

  FILE    *tau_file;
  FILE	  *vec_file;
  char    fname[80];
  int     exchange;
  double  xch_rand,tau_rand,tau_prev,tau,jump;

  int cf = check_flags(argc,argv);
  if(cf!=0){
    if(cf<0)
      printf("\n  INVALID NUMBER OF ARGUMENTS\n\n");
    else
      printf("\n  INVALID FLAG: %s\n\n",argv[cf]);
    exit(0);
  }

  double xch_sigma=xch_fwhm/C;
  double tau_sigma=tau_fwhm/C;

  printf("\n  WRITING: %.2lf GB to disk!\n\n",(double)(runs*(3*(DEC+4)))*(double)(trajectory+1)/(double)1073741824);

  /////////////////////////////////////////////////////////////////////////////
  // TOP OF MAIN LOOP

  int r;
  for(r=0;r<runs;r++){

    printf("\r  ON FILE: %d / %d",r+1,runs);
    fflush(stdout);

    if(r>0)
      seed++;
    srand(seed);

    if(print_tau){
      sprintf(fname,"rot.%d.tau",r+start);
      tau_file=fopen(fname,"w");
    }

    sprintf(fname,"rot.%d.dat",r+start);
    vec_file=fopen(fname,"w");

    if(print_head){
      fprintf(vec_file,"#\n");
      fprintf(vec_file,"# TAU MED      %lf\n",tau_median);
      fprintf(vec_file,"# TAU FWHM     %lf\n",tau_fwhm);
      fprintf(vec_file,"# TAU CORR     %lf\n",tau_corr);
      fprintf(vec_file,"#\n");
      fprintf(vec_file,"# XCH MED      %lf\n",xch_median);
      fprintf(vec_file,"# XCH FWHM     %lf\n",xch_fwhm);
      fprintf(vec_file,"# XCH CORR     %lf\n",xch_corr);
      fprintf(vec_file,"#\n");
      fprintf(vec_file,"# SEED         %d\n",seed);
      fprintf(vec_file,"#\n");
    }

    p[0]=1.0/sqrt(3.0);
    p[1]=1.0/sqrt(3.0);
    p[2]=1.0/sqrt(3.0);
    p[3]=0.0;

    fprintf(vec_file,"%+*.*lf %+*.*lf %+*.*lf\n",DEC+1,DEC,p[0],DEC+1,DEC,p[1],DEC+1,DEC,p[2]);

    tau_rand=sample_gaussian();

    int t=1;
    while(1){

      tau_prev = tau_rand;
      tau_rand = sample_gaussian();
      tau_rand = tau_corr*tau_prev + sqrt(1-tau_corr*tau_corr)*tau_rand;

      tau = pow(10.0,tau_sigma*tau_rand+tau_median);

      /***************************
      *                          *
      * PHYSICS:                 *
      *                          *
      * DIFF  = 1/(L*(L+1)*TAU)  *
      * THETA = SQRT(2*DIFF)     *
      *                          *
      ***************************/

      jump = sqrt(1.0/(3.0*tau));

      xch_rand = sample_gaussian();
      xch_rand = xch_corr*tau_rand + sqrt(1-xch_corr*xch_corr)*xch_rand;

      exchange = (int)(pow(10.0,xch_sigma*xch_rand+xch_median)+0.5);

      if(print_tau)
        fprintf(tau_file,"%lf %d\n",tau,exchange);

      int x;
      for(x=0;x<exchange;x++){

        rotate(jump*sample_rayleigh());

        fprintf(vec_file,"%+*.*lf %+*.*lf %+*.*lf\n",DEC+1,DEC,p[0],DEC+1,DEC,p[1],DEC+1,DEC,p[2]);

        if(t==trajectory)
          goto done;
        t++;

      }

    }

    done:
    fclose(vec_file);
    if(print_tau)
      fclose(tau_file);

  }

  /////////////////////////////////////////////////////////////////////////////

  printf("\n\n");

  return 0;

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////
//////  This function returns a random value
//////  sampled from a rayleigh distribution
//////  with unit width
//////

double sample_rayleigh(){
  return sqrt(-2.0*log(((double)rand()+1.0)/((double)RAND_MAX+1.0)));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////
//////  This function returns a random value
//////  sampled from a gaussian distribution
//////  with zero mean and unit variance 
//////

double sample_gaussian(){
  static int	index;
  static double	value[2];
  double	s,t,u,v;
  if(index==0){
    while(1){
      u=(double)rand()/((double)RAND_MAX/2.0)-1.0;
      v=(double)rand()/((double)RAND_MAX/2.0)-1.0;
      s=u*u+v*v;
      if(s<1.0)
        break;
    }
    t=sqrt(-2.0*log(s)/s);
    value[0]=u*t;
    value[1]=v*t;
  }
  index=(index+1)%2;
  return value[index];
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////
//////  This function rotates 4-vector {p} by 'angle'
//////

void rotate(double angle)
{

  double c[4],d[4],e,f,g;

  #if AXIS!=0

  // generate C as a random unit vector in R^3

  double h,i,j;
  e=(double)rand()/(double)RAND_MAX;
  f=sqrt(1.0-e*e);
  g=(double)(2*(rand()%2)-1);
  h=cos((double)rand()/((double)RAND_MAX/PI));
  i=sqrt(1.0-h*h);
  j=(double)(2*(rand()%2)-1);
  c[0]=e*h;
  c[1]=e*i*j;
  c[2]=f*g;

  #else

  // generate C as a random unit vector perpendicular to P

  // 1) generate unit vector A perpendicular to P

  double a[3],b[3];

  a[2]=0.0;
  a[1]=p[0];
  if(a[1]==0.0)
    a[0]=1.0;
  else{
    a[0]=-1.0*p[1];
    e=sqrt(p[0]*p[0]+p[1]*p[1]);
    a[0]/=e;
    a[1]/=e;
  }

  // 2) generate unit vector B perpendicular to both P and A

  b[0]=p[1]*a[2]-p[2]*a[1];
  b[1]=p[2]*a[0]-p[0]*a[2];
  b[2]=p[0]*a[1]-p[1]*a[0];

  // 3) generate unit vector C as a random linear combination of A and B

  e=cos((double)rand()/((double)RAND_MAX/PI));
  f=(double)(2*(rand()%2)-1);
  g=sqrt(1.0-e*e);
  c[0]=e*a[0]+f*g*b[0];
  c[1]=e*a[1]+f*g*b[1];
  c[2]=e*a[2]+f*g*b[2];

  #endif

  // turn C into quaternion conjugate

  e=sin(0.5*angle);
  c[0]*=-1.0*e;
  c[1]*=-1.0*e;
  c[2]*=-1.0*e;
  c[3]=cos(0.5*angle);

  // multiply D = PC* 

  d[0]=p[3]*c[0]+p[0]*c[3]+p[1]*c[2]-p[2]*c[1];
  d[1]=p[3]*c[1]-p[0]*c[2]+p[1]*c[3]+p[2]*c[0];
  d[2]=p[3]*c[2]+p[0]*c[1]-p[1]*c[0]+p[2]*c[3];
  d[3]=p[3]*c[3]-p[0]*c[0]-p[1]*c[1]-p[2]*c[2];

  // unconjugate C

  c[0]*=-1.0;
  c[1]*=-1.0;
  c[2]*=-1.0;

  // multiply P' = CD = CPC*

  p[0]=c[3]*d[0]+c[0]*d[3]+c[1]*d[2]-c[2]*d[1];
  p[1]=c[3]*d[1]-c[0]*d[2]+c[1]*d[3]+c[2]*d[0];
  p[2]=c[3]*d[2]+c[0]*d[1]-c[1]*d[0]+c[2]*d[3];
  p[3]=c[3]*d[3]-c[0]*d[0]-c[1]*d[1]-c[2]*d[2];

  // renormalize b/c rounding errors

  e=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  p[0]/=e;
  p[1]/=e;
  p[2]/=e;
  p[3]=0.0;

  return;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////
//////  This function checks all the parameters that
//////  the program was called with, and sets the
//////  relevant variables.
//////

int check_flags(int argc, char **argv)
{

  seed =        (int)time(0);

  runs =        1;              // number of runs
  start =       0;              // index of first run

  print_tau =   0;              // tau file print flag
  print_head =  0;              // file header print flag
  
  trajectory =  1*pow(10,5);	// trajectory length

  tau_corr =    0.0;		// correlation
  tau_fwhm =    1.0;		// of the log_10(tau_c) distribution
  tau_median =  2.0;		// of the log_10(tau_c) distribution

  xch_corr =    0.0;		// correlation
  xch_fwhm =    0.0;		// of the log_10(tau_x) distribution
  xch_median =  2.0;		// of the log_10(tau_x) distribution

  if(argc%2==0)
    return -1;

  int n;
  for(n=1;n<argc;n+=2){
    if(strcmp(argv[n],"-d")==0)
      seed=atoi(argv[n+1]);
    else if(strcmp(argv[n],"-h")==0)
      print_head=atoi(argv[n+1]);      
    else if(strcmp(argv[n],"-p")==0)
      print_tau=atoi(argv[n+1]);
    else if(strcmp(argv[n],"-r")==0)
      runs=atoi(argv[n+1]);
    else if(strcmp(argv[n],"-s")==0)
      start=atoi(argv[n+1]);
    else if(strcmp(argv[n],"-t")==0)
      trajectory=atoi(argv[n+1]);
    else if(strcmp(argv[n],"-tm")==0)
      tau_median=atof(argv[n+1]);
    else if(strcmp(argv[n],"-tf")==0)
      tau_fwhm=atof(argv[n+1]);
    else if(strcmp(argv[n],"-tc")==0)
      tau_corr=atof(argv[n+1]);
    else if(strcmp(argv[n],"-xm")==0)
      xch_median=atof(argv[n+1]);
    else if(strcmp(argv[n],"-xf")==0)
      xch_fwhm=atof(argv[n+1]);
    else if(strcmp(argv[n],"-xc")==0)
      xch_corr=atof(argv[n+1]);
    else
      return n;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////









































