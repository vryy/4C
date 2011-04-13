/*!----------------------------------------------------------------------
\file randomfield.cpp
Created on: Apr 12, 2011
\brief Class for generating samples of random fields based on spectral representation

 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/

#include "randomfield.H"
#include "mlmc.H"
#include <ctime>
#include <cstdlib>
#include <iostream>



// this was included in stopro cc
#include <cmath>
// boost currently not in use, use blitz instead
//#include <boost/random.hpp>


using namespace DRT;



/*----------------------------------------------------------------------*/
/* standard constructor */
RandomField::RandomField(unsigned int  seed,double sigma, double corr_length)
{
  // Init the necessesary stuff
  seed_ = seed;
  d_ = corr_length;
  sigma_0_ = sigma;
  pi_=M_PI;
  // The StoPro will have a period of 2*pi / Deltakappa == 2*pi*N*d / 6.
  // We want this to be >= 200.
  // ceil:= return next largest integer
  N_ = (int)ceil( 600.0 / ( pi_ * d_ ) );
  //cout << " N  " << N << endl;

   // Heuristic: PSD is of `insignificant magnitude' for
   //   abs(kappa) <= 6/d
   dkappa_ = 6.0 / (d_ * N_);
   //std::cout << "bug after line"  <<  __LINE__ << std::endl;
   //std::cout << "N ="  << N  << std::endl;
   //std::cout << "max_size of Phi  " << Phi.max_size() << std::endl;
   Phi_.reserve( N_ * N_ * N_ );
   // Boost random number generator (currently not in use
   //-------------------------------------------------------------------

   /*// This defines a random number genrator
   boost::mt19937 mt;
   // Defines random number generator for numbers between 0 and 2pi
   boost::uniform_real<double> random( 0, 2*pi_ );

   // set seed of random number generator
   // same seed produces same string of random numbers
   mt.seed(seed_);
   //try out time in seconds since janaury 1970
   //seed =(time(NULL)); */

   // Blitz random number generator
   //--------------------------------------------------------------------
   uniformclosedgen_.seed(seed_ );

   for ( int k = 0; k < N_ * N_ * N_; ++k )
     {
         //Phi_.push_back( random( mt ) );
         Phi_.push_back( uniformclosedgen_.random()*2*pi_ );

     }


}

double RandomField::EvalRandomField( double x, double y, double z)
{
  double result = 0;
    for ( int k0 = 0; k0 < N_; ++k0 )
    {
      for ( int k1 = 0; k1 < N_; ++k1 )
      {
        for ( int k2 = 0; k2 < N_; ++k2 )
        {
          result += sqrt( 2 *PowerSpectralDensity( k0 * dkappa_, k1 * dkappa_, k2 * dkappa_ ) * pow( dkappa_, 3 ) ) *
          cos( k0 * dkappa_ * x  + k1 * dkappa_ * y + k2 * dkappa_ * z +
          Phi_[k0 + k1*N_ + k2*N_*N_] );
        }
      }
    }
    return sqrt( 2 ) * result;
}

// compute power spectral density

double  RandomField::PowerSpectralDensity( double kappa_x, double kappa_y, double kappa_z )
{
  //const double pi = M_PI;
  const static double coeff = pow( sigma_0_, 2 ) * pow( d_, 3 ) / pow( 2 * sqrt( pi_ ), 3 );
  return coeff * exp( - pow( d_ * kappa_x / 2, 2 ) - pow( d_ * kappa_y / 2, 2 ) - pow( d_ * kappa_z / 2, 2 ) );
}
