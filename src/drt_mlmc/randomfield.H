/*!----------------------------------------------------------------------
\file randomfield.h
Created on: Apr 12, 2011
\brief Class for generating samples of random fields based on spectral representation

 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/

#ifdef HAVE_FFTW
#ifndef RANDOMFIELD_H_
#define RANDOMFIELD_H_
//testing
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_SerialDenseVector.h>
#include <math.h>
#include "../drt_adapter/ad_str_structure.H"

#include <cmath>
#include <ctime>

//including random number library of blitz for generation of random phase angles
#include <random/uniform.h>




class RandomField {

private:
  // some variables
  int dim_; // dimension of stochastic field
  double pi_;// = M_PI;
  double d_; //correlation length
  // set variance to 0.25 by default
  double sigma_0_; //= 0.25;

  int N_;

  double dkappa_;
  unsigned int  seed_;
  std::vector<double>  Phi_0_ ;
  std::vector<double>  Phi_1_ ;
  std::vector<double>  Phi_2_ ;
  std::vector<double>  Phi_3_ ;

  //double Phi_Array_0_ [80][80][80];
  //double Phi_Array_1_ [80][80][80];
  //double Phi_Array_2_ [80][80][80];
  //double Phi_Array_3_ [80][80][80];

  //! \brief generates uniformly distributed random numbers in interval [0;1], taken from statmech
   ranlib::UniformClosed<double> uniformclosedgen_;



public:
  // Standard constructor
  RandomField(unsigned int seed,double sigma_0, double corr_length);

  // generates new random Phase vector Phi_
  void CreateNewSample(unsigned int seed);
  // funtion to call EvalRandomfield depending on dimension
  double EvalRandomField(double x, double y, double z);
  // Evaluate RandomField at specific location
  double EvalRandomField2D( double x, double y, double z);
  double EvalRandomField3D( double x, double y, double z);

  double EvalRandomFieldCylinder( double x,double y, double z);
  //double EvalQuadrantRandomField( double x, double y, double z);

  // Calculate Power Spectral Density
  double  PowerSpectralDensity2D( double kappa_x, double kappa_y );
  double  PowerSpectralDensity3D( double kappa_x, double kappa_y, double kappa_z );

  // FFT Stuff
  double EvalRandomFieldFFT(double x, double y, double z);

}; // class RandomField

#endif /* RANDOMFIELD_H_ */
#endif // HAVE_FFTw
