/*!----------------------------------------------------------------------
\file randomfield.h
Created on: November 15, 2011
\brief Class for generating samples of gaussian and non-gaussian random fields based on spectral representation
using FFT algorithms

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
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_RCP.hpp>
#include <Epetra_SerialDenseVector.h>
#include <math.h>
//#include "../drt_adapter/ad_str_structure.H"

#include <cmath>
#include <ctime>

class GenRandomField {

private:
  //! \brief do we want/need to perform spectral matching
  bool perform_spectral_matching_;

  //! \brief which proc am I on
  int myrank_;

  //! \brief dimension of stochastic field
  int dim_;

  enum PROBABILITY_DISTRIBUTIONS
  {
    normal,
    lognormal,
    beta,
  };
  //! \brief target marginal pdf of random field
  PROBABILITY_DISTRIBUTIONS marginal_pdf_;

  enum correlation_structure_
  {
    exponential,
    gaussian,
    whatever,
  };

  //! \brief correlation length
  double d_;

  //! \briefdistribution params of target marginal distribution
  vector<double> distribution_params_;

   //! \brief sigma as calculated from pdf parameters
  double sigma_0_;

  //! \brief sigma of underlying gauss process ( changes during spectral matching)
  double sigma_ul_g_cur_it_;

  double mean_;

  // holds the discrete power spectral density of the random field
  vector<double> discrete_PSD_;

  //! \brief Simulation parameters
  //! \brief length of random field
  double periodicity_ ;

  //! \brief Number of COS terms or size of discrete_PSD
  int N_;

  //! \brief length of the FFT's
  int M_;

  //! \brief smallest frequency
  double dkappa_;

  //! \brief cutoff frequency
  double kappa_u_;

  double dx_;

  //! \brief vectors to store the random phase angles
  vector<double>  Phi_0_ ;
  vector<double>  Phi_1_ ;
  vector<double>  Phi_2_ ;
  vector<double>  Phi_3_ ;

  //! \brief/ seed for random number generator
  unsigned int  seed_;

  //! \brief parameters for evaluation i.e corners of the bounding box of the discretization
  //! \brief min bounding box
  vector<double> bb_min_;

  //! \brief max bounding box
  vector<double> bb_max_;

  //! \brief we need Pi
  double pi_;

  //! \brief  this is where we store the random field in
  double* values_;

  //! \brief define calculation method whether we use FFT or COS
  int UseFFT_;

  //! \brief Do not print to stdout
  bool reduced_output_;

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

  //! \brief simulate gaussian RandomField using FFT technique
  void SimGaussRandomFieldFFT();

  //! \brief simulate gaussian RandomField using FFT technique 3D field
  void SimGaussRandomFieldFFT3D();

  //! \brief simulate gaussian RandomField using COS technique 3D field
  double SimGaussRandomFieldCOS3D(double x, double y, double z);

  //! \brief create translation field with prescribed marginal distribution
  void TranslateToNonGaussian();

  //! \brief create translation field with prescribed marginal distribution
  void TranslateToNonGaussian( double *value);

  //! \brief helper functions for spectral matching
  double Integrate(double xmin, double xmax, double ymin, double ymax, double rho);

  //! \brief helper functions for spectral matching
  double Testfunction(double argument_x, double argument_y, double rho);

  //! \brief compute bounding box of random field from discretization
  void ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret);

}; // class GenRandomField

#endif /* RANDOMFIELD_H_ */
#endif // HAVE_FFTw
