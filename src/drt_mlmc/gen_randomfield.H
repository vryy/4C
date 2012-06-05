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
  // do we want/need to perform spectral matching
  bool perform_spectral_matching_;
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
  double sigma_0_;
  // sigma of ul gauss process ( changes during spectral matching)
  double sigma_ul_g_cur_it_;
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

  // parameters for evaluation
  vector<double> bb_min_;// min Bounding Box
  vector<double> bb_max_; // max Bounding Box
  double pi_;// = M_PI;

  // store the values
  double* values_;
  // define calculation method
  int UseFFT_;




public:
  //! \brief Standard constructor
  GenRandomField(unsigned int seed,Teuchos::RCP<DRT::Discretization> discret);

  //! \brief generates new realization of random field
  void CreateNewSample(unsigned int seed);

  //! \brief get value of random field at a specific location
  double EvalFieldAtLocation(vector<double> location, bool writetofile, bool output);

  //! \brief write random field To file
  void WriteRandomFieldToFile();

  //! \brief write PSD computed from sample to file
  void WriteSamplePSDToFile(Teuchos::RCP<Teuchos::Array <double> > sample_psd);

  //! \brief compute PSD from current sample
  void GetPSDFromSample(Teuchos::RCP<Teuchos::Array <double> > sample_psd);

  //! \brief compute PSD from current sample for 3D field
  void GetPSDFromSample3D(Teuchos::RCP<Teuchos::Array <double> > sample_psd);

  //! \brief return number of cosinus terms per dimension used for computation
  int NumberOfCosTerms() {return N_;}

  //! \brief size of array of randomfield per dimension
  int SizePerDim() {return M_;}

  //! \brief return dimension of random field (either 2 or 3)
  int Dimension() {return dim_;}




private:
  //! \brief generate new random phase angles
  void CreateNewPhaseAngles(unsigned int seed);

  //! \brief calculate discrete Power Spectral Density (PSD) based on correlation function
  void CalcDiscretePSD();

  //! \brief calculate discrete Power Spectral Density (PSD) based on correlation function for 3D fields
  void CalcDiscretePSD3D();

  //! \brief perform spectral matching of PSD if random field is non-gaussian
  void SpectralMatching();

  //! \brief perform spectral matching of PSD if random field is non-gaussian for 3D PSD's
  void SpectralMatching3D();

  //! \brief perform spectral matching of PSD if random field is non-gaussian for 3D PSD's
  void SpectralMatching3D3D();

  // simulate RandomField using FFT technique
  void SimGaussRandomFieldFFT();
  // simulate RandomField using FFT technique
  void SimGaussRandomFieldFFT3D();
  double SimGaussRandomFieldCOS3D(double x, double y, double z);
  // Create Translation field with prescribed marginal distribution
  void TranslateToNonGaussian();
  void TranslateToNonGaussian( double *value);
  // Helper Functions for spectral matching
  double Integrate(double xmin, double xmax, double ymin, double ymax, double rho);
  double Testfunction(double argument_x, double argument_y, double rho);

  // Compute Bounding Box of random field from discretization
  void ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret);


}; // class GenRandomField

#endif /* RANDOMFIELD_H_ */
#endif // HAVE_FFTw
