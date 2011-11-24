/*!----------------------------------------------------------------------
\file randomfield.h
Created on: November 15, 2011
\brief Class for generating samples of gaussian and non-gaussian random fields based on spectral representation
usign FFT algorithms

 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef HAVE_FFTW
#ifndef GENRANDOMFIELD_H_
#define GENRANDOMFIELD_H_
//testingGEN
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
#include "../drt_adapter/adapter_structure.H"

#include <cmath>
#include <ctime>

//including random number library of blitz for generation of random phase angles
#include <random/uniform.h>




class GenRandomField {

private:
  // variables of random field
  int dim_; // dimension of stochastic field
  enum PROBABILITY_DISTRIBUTIONS
  {
    normal,
    lognormal,
    beta,
  };
  PROBABILITY_DISTRIBUTIONS marginal_pdf_;
  enum correlation_structure_
  {
    exponential,
    whatever,
  };
  double d_; //correlation length
  double sigma_0_; //= 0.25;
  double mean_;
  vector<double> discrete_PSD_;
  vector<double> distribution_params_; // other distribution params

  // Simulation parameters
  double periodicity_ ;//
  int N_; // Number of cos terms
  int M_; // length of FFT's
  double dkappa_;
  double dx_;

  vector<double>  Phi_0_ ;// Vectors to store the random phase angles
  vector<double>  Phi_1_ ;
  vector<double>  Phi_2_ ;
  vector<double>  Phi_3_ ;
  unsigned int  seed_;

  // parametres for evaluation
  vector<double> bb_min_;// min Bounding Box
  vector<double> bb_max_; // max Bounding Box
  double pi_;// = M_PI;

  // store the values
  double* values_;

  //! \brief generates uniformly distributed random numbers in interval [0;1], taken from statmech
   ranlib::UniformClosed<double> uniformclosedgen_;



public:
  // Standard constructor
   GenRandomField(unsigned int seed,double sigma_0, double corr_length);
   // Compute Bounding Box of random field from discretization
   void ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret);
   double EvalFieldAtLocation(vector<double> location);

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
  //double  PowerSpectralDensity2D( double kappa_x, double kappa_y );
  //double  PowerSpectralDensity3D( double kappa_x, double kappa_y, double kappa_z );

  // FFT Stuff
  void CalcDiscretePSD();
  double EvalRandomFieldFFT(double x, double y, double z);


}; // class GenRandomField

#endif /* RANDOMFIELD_H_ */
#endif // HAVE_FFTw
#endif  // #ifdef CCADISCRET
