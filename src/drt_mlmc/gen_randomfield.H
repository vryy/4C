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
//#include "../drt_adapter/ad_str_structure.H"

#include <cmath>
#include <ctime>




class GenRandomField {

private:
  // do we need to perform spectral matching
  bool init_psd_;
  // which proc am i on
   int myrank_;
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
    gaussian,
    whatever,
  };
  double d_; //correlation length
  double sigma_0_; //= 0.25;
  // sigma of ul gauss process ( changes during spectral matching)
  double sigma_ul_g_cur_it_; //= 0.25;
  double mean_;
  vector<double> discrete_PSD_;
  vector<double> distribution_params_; // other distribution params

  // Simulation parameters
  double periodicity_ ;//
  int N_; // Number of cos terms
  int M_; // length of FFT's
  double dkappa_;
  double kappa_u_;
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




public:
  // Standard constructor
  GenRandomField(unsigned int seed,Teuchos::RCP<DRT::Discretization> discret);

  // generates new realization
  void CreateNewSample(unsigned int seed);

  // Get Value of Random Field at a specific location
  double EvalFieldAtLocation(vector<double> location, bool writetofile, bool output);

  // Write Random Field To file
  void WriteRandomFieldToFile();
  void WriteSamplePSDToFile(Teuchos::RCP<Teuchos::Array <double> > sample_psd);

  // Compute PSD from current sample
  void GetPSDFromSample(Teuchos::RCP<Teuchos::Array <double> > sample_psd);
  void GetPSDFromSample3D(Teuchos::RCP<Teuchos::Array <double> > sample_psd);

  int NumberOfCosTerms() {return N_;}

  int SizePerDim() {return M_;}

  int Dimension() {return dim_;}




private:
  // generates new random Phase vector Phi_
  void CreateNewPhaseAngles(unsigned int seed);
  // Calculate Discrete Power Spectral Density based on correlation structure
  void CalcDiscretePSD();
  void CalcDiscretePSD3D();
  // Perform Spectral Matching if RF non-gaussian
  void SpectralMatching();
  void SpectralMatching3D();
  void SpectralMatching3D3D();
  // simulate RandomField using FFT technique
  void SimGaussRandomFieldFFT();
  // simulate RandomField using FFT technique
  void SimGaussRandomFieldFFT3D();
  double TESTSimGaussRandomField3D(double x, double y, double z);
  // Create Translation field with prescribed marginal distribution
  void TranslateToNonGaussian();
  // Helper Functions for spectral matching
  double Integrate(double xmin, double xmax, double ymin, double ymax, double rho);
  double Testfunction(double argument_x, double argument_y, double rho);

  // Compute Bounding Box of random field from discretization
  void ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret);


}; // class GenRandomField

#endif /* RANDOMFIELD_H_ */
#endif // HAVE_FFTw
