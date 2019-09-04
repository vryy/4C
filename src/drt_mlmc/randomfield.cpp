/*----------------------------------------------------------------------*/
/*! \file
\brief Class for generating samples of random fields based on spectral representation

\maintainer Jonas Nitzler

\level 3
*/
/*----------------------------------------------------------------------*/

#ifdef HAVE_FFTW

#include "randomfield.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_pstream.H"

// include fftw++ stuff for multidimensional FFT
#include "fftw3.h"


#include "../drt_inpar/inpar_mlmc.H"

UQ::RandomField::RandomField(
    Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& rfp)
    : marginal_pdf_(none)
{
  // initialize some common parameters

  is_bounded_ = DRT::INPUT::IntegralValue<int>(rfp, "BOUNDED");

  rf_lower_bound_ = rfp.get<double>("LOWERBOUND");

  rf_upper_bound_ = rfp.get<double>("UPPERBOUND");

  // init largest_length
  ComputeBoundingBox(discret);

  has_spatially_variable_median_ = DRT::INPUT::IntegralValue<int>(rfp, "SPATIAL_VARIABLE_MEDIAN");
}

void UQ::RandomField::ComputeBoundingBox(Teuchos::RCP<DRT::Discretization> discret)
{
  // root bounding Box
  std::vector<double> maxrbb;

  maxrbb.push_back(-10.0e19);
  maxrbb.push_back(-10.0e19);
  maxrbb.push_back(-10.0e19);

  std::vector<double> minrbb;
  minrbb.push_back(10.0e19);
  minrbb.push_back(10.0e19);
  minrbb.push_back(10.0e19);

  bb_max_.push_back(-10.0e19);
  bb_max_.push_back(-10.0e19);
  bb_max_.push_back(-10.0e19);

  bb_min_.push_back(10.0e19);
  bb_min_.push_back(10.0e19);
  bb_min_.push_back(10.0e19);
  {
    for (int lid = 0; lid < discret->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = discret->lColNode(lid);
      // check if greater than maxrbb
      if (maxrbb[0] < node->X()[0]) maxrbb[0] = node->X()[0];
      if (maxrbb[1] < node->X()[1]) maxrbb[1] = node->X()[1];
      if (maxrbb[2] < node->X()[2]) maxrbb[2] = node->X()[2];
      // check if smaller than minrbb
      if (minrbb[0] > node->X()[0]) minrbb[0] = node->X()[0];
      if (minrbb[1] > node->X()[1]) minrbb[1] = node->X()[1];
      if (minrbb[2] > node->X()[2]) minrbb[2] = node->X()[2];
    }
  }

  discret->Comm().MaxAll(&maxrbb[0], &bb_max_[0], 3);
  discret->Comm().MinAll(&minrbb[0], &bb_min_[0], 3);

  discret->Comm().Barrier();

  // compute largest dimension of the problem
  largestlength_ = -10e19;

  for (int i = 0; i < 3; i++)
  {
    if (bb_max_[i] - bb_min_[i] > largestlength_) largestlength_ = bb_max_[i] - bb_min_[i];
  }

  if (myrank_ == 0)
  {
    IO::cout << "min " << bb_min_[0] << " " << bb_min_[1] << " " << bb_min_[2] << IO::endl;
    IO::cout << "max " << bb_max_[0] << " " << bb_max_[1] << " " << bb_max_[2] << IO::endl;
    IO::cout << "largest length " << largestlength_ << IO::endl;
  }
}


#endif  // HAVE_FFTw
