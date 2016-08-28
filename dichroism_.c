/******************************************************************************

  LINEAR DICHROISM:

    This program accepts a file listing the {x, y, z} coordinates
    of a point on the surface of a sphere as a function of time,
    and outputs the calculated (noisy) linear dichroism signal.

    It takes the input file as the only argument, and creates
    an output file of the same name with '.LD' appended.

    It is based on IDL code written by T.K.H.

******************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  COUNT_MAX	15000

extern void     camera(double*,double*,double,double,double);
extern double	sample();

int main(int argc, char *argv[]){

  FILE          *infile;
  FILE          *outfile;

  char          fname[100];

  double	steps;
  double        average;
  double        dichroism;
  double        l,r,x,y,z;

  int count;

  srand(time(0));

  if(argc==1){
    printf("\n  Please give an input file.\n\n");
    exit(0);
  }
  else{
    infile=fopen(argv[1],"r");
    if(infile==NULL){
      printf("\n  File not found.\n\n");
      exit(0);
    }
    strcpy(fname,argv[1]);
    strcat(fname,".LD");
    outfile=fopen(fname,"w");
  }
  
  ////////////////////////////////////////////////////////////////////
  //
  // STEP 1: calculate average dichroism
  //         

  steps=0.0;
  average=0.0;  
  count=0;
  while( fscanf(infile,"%lf%lf%lf",&x,&y,&z)!=EOF && count<COUNT_MAX ){
    count++;
    camera(&l,&r,x,y,z);
    steps+=1.0;
    if(l+r>0.0)
      average+=(l-r)/(l+r);
  }
  average/=steps;

  rewind(infile);

  ////////////////////////////////////////////////////////////////////
  //
  // STEP 2: output dichroism less average
  //         

  count=0;
  while( fscanf(infile,"%lf%lf%lf",&x,&y,&z)!=EOF && count<COUNT_MAX ){
    count++;
    camera(&l,&r,x,y,z);
    if(l+r>0.0)
      dichroism=(l-r)/(l+r);
    else
      dichroism=0.0;
    fprintf(outfile,"%+lf\n",average-dichroism);
  }

  ////////////////////////////////////////////////////////////////////

  fclose(infile);
  fclose(outfile);

  return 0;
}

/*****************************************************************************/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////
//////  This function calculates and stores the left and right channel
//////  camera intensities {lc, rc}, given orientation {x, y, z}.
//////
//////  The first step is to apply Fourkas.
//////
//////  The resulting intensities are of the order 10^-4, so are scaled
//////  to be of order unity.
//////
//////  (TKH did this by dividing each channel by its mean intensity,
//////  though I don't see why each channel can be scaled by a different
//////  factor.  I find that even at the shortest trajectory lengths
//////  studied, the distribution of angles is near uniform, so this
//////  funtion simply scales each channel by the intensity expected 
//////  from applying Fourkas to a uniform distribution of angles.  In
//////  doing this scaling, the excitation terms from Fourkas cancel, 
//////  as they apply to the left and right channels equally. Also,
//////  though camera intensities must be non-negative, the requirement
//////  that numerical aperatures and refractivity are between zero
//////  and one mean that taking absolute values is not needed.)
//////
//////  The scaled intensities are then multiplied by the expected number
//////  of electrons per frame and the PMT gain.
//////
//////  Next, noise is added to the intensity channels.
//////
//////  The first source is from the PMT, and scales as the sqrt of twice
//////  the intensity.  The second is due to background, and scales as a
//////  factor of the average intensity.
//////
//////  Finally, the background noise can make the intensity negative,
//////  so a constant is added to each channel, which will not affect
//////  the correlation function.  In case a rare event occurs in which
//////  the intensity is still negative, it is thresholded.
//////

void camera(double *lc, double *rc, double x, double y, double z)
{

  // CAMERA PARAMETERS

  double refractivity         = 1.000;  
  double aperture_detect      = 0.750;
  //double aperture_excite      = 0.001;
  double photons_per_frame    = 200.0;
  double photomultiplier_gain = 300.0;
  double background_amplitude = 0.300;
  double fudge_constant       = 2.000;

  double cos_theta_detect=sqrt(1.0-pow(aperture_detect/refractivity,2));

  //double cos_theta_excite=sqrt(1.0-pow(aperture_excite/refractivity,2));
          
  double a_de=(1.0/12.0)*(2.0-3.0*cos_theta_detect+pow(cos_theta_detect,3));
  double b_de=(1.0/8.0)*(cos_theta_detect-pow(cos_theta_detect,3));
  double c_de=(1.0/48.0)*(7.0-3.0*cos_theta_detect-3.0*pow(cos_theta_detect,2)-pow(cos_theta_detect,3));

  //double a_ex=(1.0/12.0)*(2.0-3.0*cos_theta_excite+pow(cos_theta_excite,3));
  //double b_ex=(1.0/8.0)*(cos_theta_excite-pow(cos_theta_excite,3));

  // GEOMETRY

  *lc  = a_de+b_de*(1-z*z)+c_de*(x*x-y*y);
  *rc  = a_de+b_de*(1-z*z)-c_de*(x*x-y*y);

  //*lc *= 2.0*a_ex+2.0*b_ex*(1-z*z);
  //*rc *= 2.0*a_ex+2.0*b_ex*(1-z*z);

  //*lc  = fabs(*lc);
  //*rc  = fabs(*rc);

  //printf("%lf   %lf   ",*lc,*rc);

  // SCALING

  *lc *= 48.0;
  *rc *= 48.0;

  *lc /= (cos_theta_detect-1.0);
  *rc /= (cos_theta_detect-1.0);

  *lc /= (pow(cos_theta_detect,2)+cos_theta_detect-8.0);
  *rc /= (pow(cos_theta_detect,2)+cos_theta_detect-8.0);

  //printf("%lf   %lf\n",*lc,*rc);

  // MEAN INTENSITY

  *lc *= photons_per_frame*photomultiplier_gain;
  *rc *= photons_per_frame*photomultiplier_gain;

  // DETECTOR NOISE

  *lc += sample()*sqrt(*lc*2.0);
  *rc += sample()*sqrt(*rc*2.0);  

  // BACKGROUND NOISE

  *lc += sample()*photons_per_frame*photomultiplier_gain*background_amplitude;
  *rc += sample()*photons_per_frame*photomultiplier_gain*background_amplitude;

  // BACKGROUND CONSTANT

  *lc += fudge_constant*photons_per_frame*photomultiplier_gain;
  *rc += fudge_constant*photons_per_frame*photomultiplier_gain;

  // THRESHOLD

  if(*lc<0.0)
    *lc=0.0;
  if(*rc<0.0)
    *rc=0.0;
              
  return;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////
//////  This function returns a random value
//////  sampled from a gaussian distribution
//////  with zero mean and unit variance. 
//////

double sample()
{
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










