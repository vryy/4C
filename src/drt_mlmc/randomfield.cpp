/*!----------------------------------------------------------------------
\file randomfield_fft.cpp
Created on: November, 2011
\brief Class for generating samples of random fields based on spectral representation
using FFT algorithms

 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/
#ifdef HAVE_FFTW

#include "randomfield.H"
#include "../drt_lib/drt_globalproblem.H"
#include "mlmc.H"
#include <ctime>
#include <cstdlib>
#include <iostream>

#include <complex>



// this was included in stopro cc
#include <cmath>
// boost currently not in use, use blitz instead
#include <boost/random.hpp>
// include fftw++ tsuff for multidimensional FFT
#include"fftw3.h"
//include "/home/biehler/software/fft_stuff/fftw++-1.09/fftw++.h"
//#include "/home/biehler/software/fft_stuff/fftw++-1.09/Array.h"
//using namespace Array;
//%using namespace fftwpp;
using namespace DRT;



/*----------------------------------------------------------------------*/
/* standard constructor */
RandomField::RandomField(unsigned int  seed,double sigma, double corr_length)
{

  // Init the necessesary stuff
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  // Dimension
  dim_ = mlmcp.get<int>("RANDOM_FIELD_DIMENSION");
  if(dim_!=3&&dim_!=2)
      dserror("Dimension of random field must be 2 or 3, fix your input file");
  N_= mlmcp.get<int>("NUM_COS_TERMS");
  seed_ = seed;
  d_ = corr_length;
  sigma_0_ = sigma;
  pi_=M_PI;
  // The StoPro will have a period of 2*pi / Deltakappa == 2*pi*N*d / 6.
  // We want this to be >= 200.
  // ceil:= return next largest integer
  //N_ = (int)ceil( 600.0 / ( pi_ * d_ ) );
// Heuristic: PSD is of `insignificant magnitude' for
   //   abs(kappa) <= 6/d

  dkappa_ = 2*pi_/(mlmcp.get<double>("PERIODICITY")*2);

   // Boost random number generator (currently not in use
   //-------------------------------------------------------------------
  // This defines a random number genrator
   boost::mt19937 mt;
   // Defines random number generator for numbers between 0 and 2pi
   boost::uniform_real<double> random( 0, 2*pi_ );
   mt.seed(seed_);


   //ofstream File("OutputPhi.csv");

   switch(dim_){
   case 3:
     Phi_0_.reserve( N_ * N_ * N_ );
     Phi_1_.reserve( N_ * N_ * N_ );
     Phi_2_.reserve( N_ * N_ * N_ );
     Phi_3_.reserve( N_ * N_ * N_ );

     for ( int k5 = 0; k5 < N_ * N_ * N_; ++k5 )
     {
       Phi_0_.push_back( random( mt ) );
       Phi_1_.push_back( random( mt ) );
       Phi_2_.push_back( random( mt ) );
       Phi_3_.push_back( random( mt ) );
     }
     break;
   case 2:
        Phi_0_.reserve( N_ * N_ );
        Phi_1_.reserve( N_ * N_ );
        for ( int k5 = 0; k5 < N_ * N_; ++k5 )
        {
          Phi_0_.push_back( random( mt ) );
          Phi_1_.push_back( random( mt ) );
        }
        break;
   default:
     dserror("Dimension of random field must be 2 or 3, fix your input file");
     break;
   }
}
void RandomField::CreateNewSample(unsigned int seed)
{

  // This defines a random number genrator
  boost::mt19937 mt;
  // Defines random number generator for numbers between 0 and 2pi
  boost::uniform_real<double> random( 0, 2*pi_ );

  // set seed of random number generator
  // same seed produces same string of random numbers
  mt.seed(seed);
  //try out time in seconds since janaury
  switch (dim_){

  case 3:
    Phi_0_.clear();
    Phi_1_.clear();
    Phi_2_.clear();
    Phi_3_.clear();

    for ( int k5 = 0; k5 < N_ * N_ * N_; ++k5 )
      {
        Phi_0_.push_back( random( mt ) );
        Phi_1_.push_back( random( mt ) );
        Phi_2_.push_back( random( mt ) );
        Phi_3_.push_back( random( mt ) );
        //.push_back( random( mt ) );
      }
    break;
  case 2:
    Phi_0_.clear();
    Phi_1_.clear();
    for ( int k5 = 0; k5 < N_ * N_ ; ++k5 )
         {
           Phi_0_.push_back( random( mt ) );
           Phi_1_.push_back( random( mt ) );

         }
    break;
  default:
    dserror("Dimension of random field must be 2 or 3, fix your input file");
    break;
  }

}

// Wrapper to decide wether we are 2D or 3D
double RandomField::EvalRandomField(double x, double y, double z)
{
  double result=0;
  switch (dim_)
  {
  case 3:
    result=EvalRandomField3D(x,y,z);
    break;
  case 2:
    result=EvalRandomField2D(x,y,z);
    break;
  }
  return result;
}


double RandomField::EvalRandomField2D( double x, double y, double z)

{
  double result = 0;
    for ( int k0 = 0; k0 < N_; ++k0 )
    {
      for ( int k1 = 0; k1 < N_; ++k1 )
      {
       /* result += sqrt( 2*PowerSpectralDensity2D( (k0+1) * dkappa_, (k1+1) * dkappa_) * pow( dkappa_, 2 ) ) *
                  (cos( (k0+1) * dkappa_ * x  + (k1+1) * dkappa_ * y  + Phi_0_[k0 + k1*N_])+
                      cos( (k0+1) * dkappa_ * x  - (k1+1) * dkappa_ * y  + Phi_1_[k0 + k1*N_]) );*/
        // test forsumm from k=0
        result += sqrt( 2*PowerSpectralDensity2D( (k0) * dkappa_, (k1) * dkappa_) * pow( dkappa_, 2 ) ) *
                          (cos( (k0) * dkappa_ * x  + (k1) * dkappa_ * y  + Phi_0_[k0 + k1*N_])+
                              cos( (k0) * dkappa_ * x  - (k1) * dkappa_ * y  + Phi_1_[k0 + k1*N_]) );
      }
    }

    return sqrt( 2 ) * result;
}


double RandomField::EvalRandomField3D( double x, double y, double z)

{
  //ofstream File("OutputPhiIndex.csv");
  double result = 0;
    for ( int k0 = 0; k0 < N_; ++k0 )
    {
      for ( int k1 = 0; k1 < N_; ++k1 )
      {
        for ( int k2 = 0; k2 < N_; ++k2 )
        {
          // This tiny little formular is from Shinozuka1996, page 46 Simulation of quadrant fields
          result += sqrt( 2*PowerSpectralDensity3D( (k0+1) * dkappa_, (k1+1) * dkappa_, (k2+1) * dkappa_ ) * pow( dkappa_, 3 ) ) *
                    (cos( (k0+1) * dkappa_ * x  + (k1+1) * dkappa_ * y + (k2+1) * dkappa_ * z + Phi_0_[k0 + k1*N_ + k2*N_*N_])+
                        cos( (k0+1) * dkappa_ * x  + (k1+1) * dkappa_ * y - (k2+1) * dkappa_ * z + Phi_1_[k0 + k1*N_ + k2*N_*N_])+
                        cos( (k0+1) * dkappa_ * x  - (k1+1) * dkappa_ * y - (k2+1) * dkappa_ * z + Phi_2_[k0 + k1*N_ + k2*N_*N_])+
                        cos( (k0+1) * dkappa_ * x  - (k1+1) * dkappa_ * y + (k2+1) * dkappa_ * z + Phi_3_ [k0 + k1*N_ + k2*N_*N_] )
                        );
        }
      }
    }
    return sqrt( 2 ) * result;
}

double RandomField::EvalRandomFieldCylinder( double x,double y, double z)
{
  // Instead of calculating a true 3D random field calculate
  // s = r *theta , the curvelength on the circle

  // define centerline/centerpoint of circel
  double x_center = 0;
  double y_center = 0;

  // calc r vector
  vector <double> r;

  r.push_back(x-x_center);
  r.push_back(y-y_center);

  //calc angel theta
  double theta =acos(r[0]/sqrt(pow(r[0],2)+pow(r[1],2)));
  double s = theta*sqrt(pow(r[0],2)+pow(r[1],2));
  double result = EvalRandomField2D(s, z, z);
  return result;


}


// compute power spectral density

double  RandomField::PowerSpectralDensity3D( double kappa_x, double kappa_y, double kappa_z )
{
  const static double coeff = pow( sigma_0_, 2 ) * pow( d_, 3 ) / pow( (2 * sqrt( pi_ )), 3 );
  return coeff * exp( - pow( d_ * kappa_x / 2, 2 ) - pow( d_ * kappa_y / 2, 2 ) - pow( d_ * kappa_z / 2, 2 ) );
}

double  RandomField::PowerSpectralDensity2D( double kappa_x, double kappa_y)
{
  const static double coeff = pow( sigma_0_, 2 ) * pow( d_, 2 ) / pow( (2 * sqrt( pi_ )), 2 );
  return coeff * exp( - pow( d_ * kappa_x / 2, 2 ) - pow( d_ * kappa_y / 2, 2 ) );
}

// HERE comes the experimental FFT Stuff
double RandomField::EvalRandomFieldFFT(double x, double y, double z)
{
// Lets see if we can speed up things with FFTW
  int M =N_*2; // Define number of points
  // double for loops to compute coefficients
  double A; // store some stuff
  // store coefficients


  complex<double>* b1;
  complex<double>* b2;
  b1 = new complex<double>[M*M];
  b2 = new complex<double>[M*M];
  // define complex i
  complex<double> i_comp (0,1);


  for (int j=0;j<M;j++)
  {
    for (int k=0;k<M;k++)
    {
      //A=sqrt(2*(pow(sigma_0_,2)*d_/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
      // sort entries ro w major style
      // set first elements to zero
      if(k==0||j==0||j>(N_-2)||k>(N_-2))
      {
        b1[k+M*j]=0.0;
        b2[k+M*j]=0.0;
      }
      else
      {
        A=sqrt(2*(pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/2,2))-(pow(d_*k*dkappa_/2,2))))*(pow(dkappa_,2)));
        real(b1[k+M*j])=A*sqrt(2)*cos(Phi_1_[k+M*j]);
        imag(b1[k+M*j])= A*sqrt(2)*sin(Phi_1_[k+M*j]);
        real(b2[k+M*j])= A*sqrt(2)*cos(Phi_2_[k+M*j]);
        imag(b2[k+M*j])= A*sqrt(2)*sin(Phi_2_[k+M*j]);
      }
     }
  }



  int rank = 1; /* not 2: we are computing 1d transforms */
  //int n[] = {1024}; /* 1d transforms of length 10 */
  int N_fftw = M;
  int howmany = M; // same here
  int idist = M;
  int    odist = M;
  int istride =1;
  int ostride = 1; /* distance between two elements in the same column */
  //int *inembed = n;
  //int *onembed = n;

  // allocate output arrays
  complex<double>* d1;
  complex<double>* d2;
  d1 = new complex<double>[M*M];
  d2 = new complex<double>[M*M];
  fftw_plan ifft_of_rows;
  fftw_plan ifft_of_rows2;
  fftw_plan ifft_of_collums;

  ifft_of_rows = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(b1)),
      NULL,istride, idist,(reinterpret_cast<fftw_complex*>(d1)),
                  NULL,
                  ostride,
                  odist,
                  FFTW_BACKWARD,
                  FFTW_ESTIMATE);
  ifft_of_rows2 = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(b2)),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(d2)),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);
  //ifft_of_rows = fftw_plan_dft_1d(1024, (reinterpret_cast<fftw_complex*>(b1)), (reinterpret_cast<fftw_complex*>(b1)), FFTW_FORWARD, FFTW_ESTIMATE);

  istride =M;
  ostride=M;
  idist=1;
  odist=1;
  ifft_of_collums = fftw_plan_many_dft(rank, &N_fftw,howmany,(reinterpret_cast<fftw_complex*>(d1)),
       NULL,istride, idist,(reinterpret_cast<fftw_complex*>(d2)),
                   NULL,
                   ostride,
                   odist,
                   FFTW_BACKWARD,
                   FFTW_ESTIMATE);

  // just execute for the hack of it
  // transpose d1


  fftw_execute(ifft_of_rows);
  //for (int k=0;k<M*9;k++)
    // {
     // cout<< "d1" << d1[k] << endl;
    // }
     //dserror("stop right here");

  fftw_execute(ifft_of_rows2);

complex<double> scaling (M,M);
  // transpose d1
   for (int k=0;k<M*M;k++)
   {
     d2[k]=conj(d2[k]);
     d1[k]=d1[k]+d2[k];
   }


   fftw_execute(ifft_of_collums);


   // get coordinates
   int index_x=int(floor(x));
   int index_y=int(floor(y));
   //cout << "index_x :" << index_x <<endl;
   //cout << "index_y :" << index_y <<endl;
   double youngs = real(d2[index_x+index_y*128]);
   //cout <<" youngs "<< youngs << endl;
   // free memory
   delete b1;
   delete b2;
   delete d1;
   delete d2;
   fftw_destroy_plan(ifft_of_rows);
   fftw_destroy_plan(ifft_of_collums);


return youngs;


}
#endif
