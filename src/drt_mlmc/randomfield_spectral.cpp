/*----------------------------------------------------------------------*/
/*! \file
\brief Class for generating samples of gaussian and non gaussian random fields based on spectral
representation


\level 2
*/
/*----------------------------------------------------------------------*/

#ifdef HAVE_FFTW

// disable boost long double functions because they do not work together with
// valgrind. valgrind emulates long double as double and therefore some
// routines inside boost do not converge to the expected precision, crashing
// the program at the global initialization stage (before main()). we never
// use <long double> distributions anyway, so no harm is caused by disabling
// these methods.
#define BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS


#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "randomfield.H"
#include "randomfield_spectral.H"
#include "mlmc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_pstream.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <complex>
#include <cmath>
// For colored std::couts
#include "../drt_lib/drt_colors.H"
#include <boost/random.hpp>

// include fftw++ stuff for multidimensional FFT
#include "fftw3.h"

#include <fstream>

#include <boost/math/distributions/beta.hpp>       // for beta_distribution.
#include <boost/math/distributions/normal.hpp>     // for normal_distribution.
#include <boost/math/distributions/lognormal.hpp>  // for lognormal _distribution.
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/accumulators.hpp>


#include "../drt_inpar/inpar_mlmc.H"

using boost::accumulators::mean;
using boost::accumulators::stats;
using boost::accumulators::variance;
using boost::math::beta_distribution;
using boost::math::lognormal_distribution;
using boost::math::normal_distribution;

/*----------------------------------------------------------------------*/
/* standard constructor */
UQ::RandomFieldSpectral::RandomFieldSpectral(
    unsigned int seed, Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& rfp)
    : RandomField(discret, rfp)
{
  myrank_ = discret->Comm().MyPID();

  // spatial dimension of random field  only 2 adn 3 supported
  dim_ = rfp.get<int>("RANDOM_FIELD_DIMENSION");
  if (dim_ != 3 && dim_ != 2)
    dserror("Dimension of random field must be 2 or 3, fix your input file");

  // do we want to perform spectral matching with PSd
  perform_spectral_matching_ = DRT::INPUT::IntegralValue<int>(rfp, "SPECTRAL_MATCHING");
  N_ = rfp.get<int>("NUM_COS_TERMS");
  kappa_u_ = rfp.get<double>("KAPPA_U");
  seed_ = seed;
  d_ = rfp.get<double>("CORRLENGTH");
  sigma_0_ = rfp.get<double>("SIGMA");
  sigma_ul_g_cur_it_ = 0.0;
  pi_ = M_PI;
  M_ = rfp.get<int>("SIZE_PER_DIM");
  dkappa_ = kappa_u_ / N_;
  periodicity_ = 2. * pi_ / dkappa_;
  dx_ = periodicity_ / M_;

  if (periodicity_ < 1.1 * largestlength_)
    dserror("Periodic length of random field is to small for your problem ");

  if (periodicity_ > 2.0 * largestlength_)
  {
    IO::cout << "periodicity_ " << periodicity_ << IO::endl;
    IO::cout << " WARNING Periodic length of random field is very large compared to your "
                "discretization, you should consider a smaller value"
             << IO::endl;
  }
  det_value_paracont_ = rfp.get<double>("CONTBLENDVALUE");

  // distribution parameters of non gaussian pdf
  distribution_params_.push_back(rfp.get<double>("NONGAUSSPARAM1"));
  distribution_params_.push_back(rfp.get<double>("NONGAUSSPARAM2"));

  // Get correlation structure so far we can only use gaussian
  INPAR::MLMC::CorrStruct cstruct =
      DRT::INPUT::IntegralValue<INPAR::MLMC::CorrStruct>(rfp, "CORRSTRUCT");
  switch (cstruct)
  {
    case INPAR::MLMC::corr_gaussian:
      // nothing to do here
      break;
    default:
      dserror("Unknown Correlation structure");
      break;
  }
  // do we want to truncate the field


  //
  double upper_bound;

  // compute parameters for non gaussian pdf
  INPAR::MLMC::MarginalPdf mpdf =
      DRT::INPUT::IntegralValue<INPAR::MLMC::MarginalPdf>(rfp, "MARGINALPDF");
  switch (mpdf)
  {
    case INPAR::MLMC::pdf_gaussian:
      marginal_pdf_ = normal;
      // Hack
      IO::cout << "remove a line here" << IO::endl;
      sigma_0_ = sqrt((exp(pow(distribution_params_[1], 2)) - 1) *
                      exp(2 * distribution_params_[0] + pow(distribution_params_[1], 2)));
      break;
    case INPAR::MLMC::pdf_beta:
      marginal_pdf_ = beta;
      // compute bounds of distribution
      // using mu_b = 0 and sigma_b = 1 the lower and upper bounds can be computed according to
      // Yamazaki1988
      // lower bound
      distribution_params_.push_back(
          -1.0 *
          sqrt(distribution_params_[0] * (distribution_params_[0] + distribution_params_[1] + 1) /
               distribution_params_[1]));
      // upper bound
      upper_bound = (-1.0 * sqrt(distribution_params_[1] *
                                 (distribution_params_[0] + distribution_params_[1] + 1) /
                                 distribution_params_[0]));
      // I need abs(lB)+abs(uB) hence
      distribution_params_.push_back(abs(upper_bound) + abs(distribution_params_[2]));
      IO::cout << "Distribution parameter of beta distribution " << distribution_params_[0] << " "
               << distribution_params_[1] << " " << distribution_params_[2] << " "
               << distribution_params_[3] << IO::endl;
      break;
    case INPAR::MLMC::pdf_lognormal:
      marginal_pdf_ = lognormal;
      // Calculate mean of lognormal distribution based mu_N and sigma_N
      distribution_params_.push_back(
          exp(distribution_params_[0] + 0.5 * pow(distribution_params_[1], 2)));
      // also calc variance and sigma
      //(exp(s^2 )-1) * exp(2m + s^2)
      sigma_0_ = sqrt((exp(pow(distribution_params_[1], 2)) - 1) *
                      exp(2 * distribution_params_[0] + pow(distribution_params_[1], 2)));

      if (myrank_ == 0)
      {
        IO::cout << "sigma_0 " << sigma_0_ << IO::endl;
        IO::cout << "Distribution parameter of lognormal distribution " << distribution_params_[0]
                 << " " << distribution_params_[1] << " " << distribution_params_[2] << IO::endl;
      }
      break;
    default:
      dserror("Unknown Marginal pdf");
      break;
  }

  // Get calculation method

  INPAR::MLMC::CalcMethod calcm =
      DRT::INPUT::IntegralValue<INPAR::MLMC::CalcMethod>(rfp, "CALC_METHOD");
  switch (calcm)
  {
    case INPAR::MLMC::calc_m_fft:
      UseFFT_ = 1;
      break;
    case INPAR::MLMC::calc_m_cos:
      UseFFT_ = 0;
      break;
    default:
      dserror("Unknown Calculation Method for RF choose fft or cos");
      break;
  }

  // create Multidimesional array to store the values
  double mydim = dim_;
  // transform needed because pow does not like to get two ints
  int size_of_field = int(pow(M_, mydim));
  values_ = new double[size_of_field];
  // The random field will have a period of 2*pi / Deltakappa == 2*pi*N*d / 6.

  if (myrank_ == 0)
  {
    IO::cout << "Random Field Parameters " << IO::endl;
    IO::cout << "Periodicity L: " << periodicity_ << IO::endl;
    IO::cout << "M: " << M_ << IO::endl;
    IO::cout << "N: " << N_ << IO::endl;
    IO::cout << "kappa_u: " << kappa_u_ << IO::endl;
    IO::cout << "dkappa: " << dkappa_ << IO::endl;
    IO::cout << "dx " << dx_ << IO::endl;
  }

  // reserve some space to store the random phase angles in
  switch (dim_)
  {
    case 3:
      Phi_0_.reserve(N_ * N_ * N_);
      Phi_1_.reserve(N_ * N_ * N_);
      Phi_2_.reserve(N_ * N_ * N_);
      Phi_3_.reserve(N_ * N_ * N_);
      discrete_PSD_.reserve(int(pow(N_, 3.0)));
      break;
    case 2:
      Phi_0_.reserve(N_ * N_);
      Phi_1_.reserve(N_ * N_);
      discrete_PSD_.reserve(int(pow(N_, 2.0)));
      break;
    default:
      dserror("Dimension of random field must be 2 or 3, fix your input file");
      break;
  }
  // ComputeBoundingBox(discret);
  CreateNewPhaseAngles(seed_);

  switch (dim_)
  {
    case 3:
      CalcDiscretePSD3D();
      if (UseFFT_) SimGaussRandomFieldFFT3D();
      break;
    case 2:
      CalcDiscretePSD();
      if (UseFFT_) SimGaussRandomFieldFFT();
      break;
    default:
      dserror("Dimension of random field must be 2 or 3, fix your input file");
      break;
  }
  // WriteRandomFieldToFile();
  if (UseFFT_) TranslateToNonGaussian();
}
void UQ::RandomFieldSpectral::CreateNewSample(unsigned int seed)
{
  // check wether we use fft
  CreateNewPhaseAngles(seed);
  if (UseFFT_)
  {
    switch (dim_)
    {
      case 3:
        SimGaussRandomFieldFFT3D();
        break;
      case 2:
        SimGaussRandomFieldFFT();
        ;
        break;
      default:
        dserror("Dimension of random field must be 2 or 3, fix your input file");
        break;
    }
    TranslateToNonGaussian();
  }
}


void UQ::RandomFieldSpectral::CreateNewPhaseAngles(unsigned int seed)
{
  // This defines a random number genrator
  boost::mt19937 mt;
  // Defines random number generator for numbers between 0 and 2pi
  boost::uniform_real<double> random(0, 2 * pi_);

  // set seed of random number generator
  // same seed produces same std::string of random numbers // on the same platform :-)
  mt.seed(seed);
  switch (dim_)
  {
    case 3:
      Phi_0_.clear();
      Phi_1_.clear();
      Phi_2_.clear();
      Phi_3_.clear();

      for (int k5 = 0; k5 < N_ * N_ * N_; ++k5)
      {
        Phi_0_.push_back(random(mt));
        Phi_1_.push_back(random(mt));
        Phi_2_.push_back(random(mt));
        Phi_3_.push_back(random(mt));
      }
      break;
    case 2:
      Phi_0_.clear();
      Phi_1_.clear();
      for (int k5 = 0; k5 < N_ * N_; ++k5)
      {
        Phi_0_.push_back(random(mt));
        Phi_1_.push_back(random(mt));
      }
      break;
    default:
      dserror("Dimension of random field must be 2 or 3, fix your input file");
      break;
  }
  if (myrank_ == 0)
  {
    std::ofstream File;
    File.open("Phi.txt", std::ios::out);
    int size = int(Phi_0_.size());
    for (int i = 0; i < size; i++)
    {
      File << std::setprecision(15) << Phi_0_[i] << std::endl;
    }
    File.close();
  }
}

// compute power spectral density
void UQ::RandomFieldSpectral::CalcDiscretePSD()
{
  // just compute PSD
  // IO::cout<< "sigma_0_" << sigma_0_ << IO::endl;
  for (int j = 0; j < N_; j++)
  {
    for (int k = 0; k < N_; k++)
    {
      if (k == 0 || j == 0)
      {
        discrete_PSD_.push_back(
            0.5 * (pow(sigma_0_, 2) * pow(d_, 2) / (4 * pi_) *
                      exp(-(pow(d_ * j * dkappa_ / 2, 2)) - (pow(d_ * k * dkappa_ / 2, 2)))));
      }
      else
      {
        // discrete_PSD_.push_back((pow(sigma_0_,2)*pow(d_,2)/(4*pi_)*exp(-(pow(d_*j*dkappa_/(2*sqrt(pi_)),2))-(pow(d_*k*dkappa_/(2*sqrt(pi_)),2)))));
        discrete_PSD_.push_back(
            (pow(sigma_0_, 2) * pow(d_, 2) / (4 * pi_) *
                exp(-(pow(d_ * j * dkappa_ / 2, 2)) - (pow(d_ * k * dkappa_ / 2, 2)))));
      }
    }
  }
  // Write to file
  /* if (myrank_ == 0)
      {
        std::ofstream File;
        File.open("DiscretePSDnewstyle.txt",std::ios::out);
        int size = int (pow(N_,2.0));
        for(int i=0;i<size;i++)
        {
          File << discrete_PSD_[i]<< std::endl;
        }
        File.close();
      }*/

  if (marginal_pdf_ != normal)
  {
    if (perform_spectral_matching_)
    {
      SpectralMatching();
    }
    else
    {
      IO::cout << RED_LIGHT " WARNING NO SPECTRAL MATCHING BUT PDF NONGAUSS " END_COLOR << IO::endl;
      sigma_ul_g_cur_it_ = sigma_0_;
    }  // compute underlying gaussian distribution based on shields2011
  }
  else
  {
    if (myrank_ == 0) IO::cout << " Nothing to do marginal pdf gaussian " << IO::endl;
  }
}

void UQ::RandomFieldSpectral::CalcDiscretePSD3D()
{
  // just compute PSD
  for (int j = 0; j < N_; j++)
  {
    for (int k = 0; k < N_; k++)
    {
      for (int l = 0; l < N_; l++)
      {
        if (k == 0 || j == 0 || l == 0)
        {
          discrete_PSD_.push_back(
              0.5 * (pow(sigma_0_, 2) * pow(d_, 3) / (pow((2 * sqrt(pi_)), 3)) *
                        exp(-(pow(d_ * j * dkappa_ / 2, 2)) - (pow(d_ * k * dkappa_ / 2, 2)) -
                            (pow(d_ * l * dkappa_ / 2, 2)))));
        }
        else
        {
          discrete_PSD_.push_back(
              (pow(sigma_0_, 2) * pow(d_, 3) / (pow((2 * sqrt(pi_)), 3)) *
                  exp(-(pow(d_ * j * dkappa_ / 2, 2)) - (pow(d_ * k * dkappa_ / 2, 2)) -
                      (pow(d_ * l * dkappa_ / 2, 2)))));
        }
      }
    }
  }
  //  // Write to file
  /*  if (myrank_ == 0)
    {
      std::ofstream File;
      File.open("DiscretePSDnewstyle.txt",std::ios::out);
      int size = int (pow(N_,3.0));
      for(int i=0;i<size;i++)
      {
        File << discrete_PSD_[i]<< std::endl;
      }
      File.close();
    }*/


  if (marginal_pdf_ != normal)
  {
    // compute underlying gaussian distribution based on shields2011
    // IO::cout<< "NO SPECTRAL MATCHING"<< IO::endl;
    if (perform_spectral_matching_)
    {
      SpectralMatching3D();
      // SpectralMatching3D3D();
    }
    else
    {
      IO::cout << RED_LIGHT " WARNING NO SPECTRAL MATCHING BUT PDF NONGAUSS " END_COLOR << IO::endl;
      sigma_ul_g_cur_it_ = sigma_0_;
    }
  }
  else
  {
    if (myrank_ == 0) IO::cout << " Nothing to do marginal pdf gaussian " << IO::endl;
  }
}


void UQ::RandomFieldSpectral::SimGaussRandomFieldFFT()
{
  double A;  // store some stuff

  Teuchos::RCP<Teuchos::Array<std::complex<double>>> b1 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> b2 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_));

  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d1 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d2 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_));


  std::complex<double> i_comp(0, 1);
  for (int j = 0; j < M_; j++)
  {
    for (int k = 0; k < M_; k++)
    {
      // sort entries row major style set first elements to zero
      if (k == 0 || j == 0 || j > (N_ - 2) || k > (N_ - 2))
      {
        (*b1)[k + M_ * j] = 0.0;
        (*b2)[k + M_ * j] = 0.0;
      }
      else
      {
        A = sqrt(2 * (discrete_PSD_[k + j * N_] * (pow(dkappa_, 2))));
        ((*b1)[k + M_ * j]).real(A * sqrt(2) * cos(Phi_0_[k + N_ * j]));
        ((*b1)[k + M_ * j]).imag(A * sqrt(2) * sin(Phi_0_[k + N_ * j]));
        ((*b2)[k + M_ * j]).real(A * sqrt(2) * cos(Phi_1_[k + N_ * j]));
        ((*b2)[k + M_ * j]).imag(A * sqrt(2) * sin(Phi_1_[k + N_ * j]));
      }
    }
  }


  int rank = 1; /* not 2: we are computing 1d transforms 1d transforms of length M_ */
  int N_fftw = M_;
  int howmany = M_;  // same here
  int idist = M_;    //  the distance in memory between the first element of the first array and the
                     //  first element of the second array */
  int odist = M_;
  int istride = 1;
  int ostride = 1; /* distance between two elements in the same row/column/rank */
  // int *inembed = n;
  // int *onembed = n;


  fftw_plan ifft_of_rows;
  fftw_plan ifft_of_rows2;
  fftw_plan ifft_of_collums;
  // Why we have to do this &((*b1)[0]))
  // the fftw requires the adress of the first element of a fftw_complex array. We can  get to a
  // normal array by casting. Since we do not want to use standard arrays with new double[size] but
  // Teuchos::RCP::Teuchos::Arrays instead we have to do (*b1)[0] to get the first element of the
  // array and than &(..) to get the adress

  ifft_of_rows =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*b1)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d1)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  ifft_of_rows2 =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*b2)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d2)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);

  istride = M_;
  ostride = M_;
  idist = 1;
  odist = 1;

  ifft_of_collums =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d2)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(ifft_of_rows);
  fftw_execute(ifft_of_rows2);
  // transpose d1
  for (int k = 0; k < M_ * M_; k++)
  {
    (*d2)[k] = conj((*d2)[k]);
    (*d1)[k] = (*d1)[k] + (*d2)[k];
  }
  fftw_execute(ifft_of_collums);

  for (int i = 0; i < M_ * M_; i++)
  {
    values_[i] = ((*d2)[i]).real();
  }
  // free memory
  fftw_destroy_plan(ifft_of_rows);
  fftw_destroy_plan(ifft_of_rows2);
  fftw_destroy_plan(ifft_of_collums);
}
void UQ::RandomFieldSpectral::SimGaussRandomFieldFFT3D()
{
  double A;  // store some stuff
  // store coefficients
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> b1 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> b2 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> b3 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> b4 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d1 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d2 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d3 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d4 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d5 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d6 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d7 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));

  std::complex<double> i_comp(0, 1);
  for (int j = 0; j < M_; j++)
  {
    for (int k = 0; k < M_; k++)
    {
      for (int l = 0; l < M_; l++)
      {
        // sort entries row major style set first elements to zero
        if (k == 0 || j == 0 || l == 0 || j > (N_ - 2) || k > (N_ - 2) || l > (N_ - 2))
        {
          (*b1)[l + M_ * (k + M_ * j)] = 0.0;
          (*b2)[l + M_ * (k + M_ * j)] = 0.0;
          (*b3)[l + M_ * (k + M_ * j)] = 0.0;
          (*b4)[l + M_ * (k + M_ * j)] = 0.0;
        }
        else
        {
          A = sqrt(2 * (discrete_PSD_[l + N_ * (k + N_ * j)] * (pow(dkappa_, 3))));
          ((*b1)[l + M_ * (k + M_ * j)]).real(A * sqrt(2) * cos(Phi_0_[l + N_ * (k + N_ * j)]));
          ((*b1)[l + M_ * (k + M_ * j)]).imag(A * sqrt(2) * sin(Phi_0_[l + N_ * (k + N_ * j)]));
          ((*b2)[l + M_ * (k + M_ * j)]).real(A * sqrt(2) * cos(Phi_1_[l + N_ * (k + N_ * j)]));
          ((*b2)[l + M_ * (k + M_ * j)]).imag(A * sqrt(2) * sin(Phi_1_[l + N_ * (k + N_ * j)]));
          ((*b3)[l + M_ * (k + M_ * j)]).real(A * sqrt(2) * cos(Phi_2_[l + N_ * (k + N_ * j)]));
          ((*b3)[l + M_ * (k + M_ * j)]).imag(A * sqrt(2) * sin(Phi_2_[l + N_ * (k + N_ * j)]));
          ((*b4)[l + M_ * (k + M_ * j)]).real(A * sqrt(2) * cos(Phi_3_[l + N_ * (k + N_ * j)]));
          ((*b4)[l + M_ * (k + M_ * j)]).imag(A * sqrt(2) * sin(Phi_3_[l + N_ * (k + N_ * j)]));
        }
      }
    }
  }
  // set up FFTW Plans
  int rank = 1; /* not 2: we are computing 1d transforms  1d transforms of length M_*/
  int N_fftw = M_;
  int howmany = M_ * M_;  // same here
  int idist = M_;  //  the distance in memory between the first element of the first array and the
                   //  first element of the second array */
  int odist = M_;
  int istride = 1;
  int ostride = 1; /* distance between two consecutive elements */
  // four plans for the rows
  fftw_plan ifft_of_rows1;
  fftw_plan ifft_of_rows2;
  fftw_plan ifft_of_rows3;
  fftw_plan ifft_of_rows4;
  // two for the collumns
  fftw_plan ifft_of_collums1;
  fftw_plan ifft_of_collums2;
  // and one for the rank ( third dim)
  fftw_plan ifft_of_rank;

  ifft_of_rows1 =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*b1)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d1)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  ifft_of_rows2 =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*b2)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d2)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  ifft_of_rows3 =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*b3)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d3)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  ifft_of_rows4 =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*b4)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d4)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  // start here we need to transpose the arrays first
  // stride and dist are the same for collumns because arrays are transposed before FFTs
  istride = 1;
  ostride = 1;
  idist = M_;
  odist = M_;
  ifft_of_collums1 =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d5)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  ifft_of_collums2 =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*d3)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d6)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  // and now the rank
  istride = M_ * M_;
  ostride = M_ * M_;
  idist = 1;
  odist = 1;
  ifft_of_rank =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*d5)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d7)[0]))), NULL, ostride,
          odist, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(ifft_of_rows1);
  fftw_execute(ifft_of_rows2);
  fftw_execute(ifft_of_rows3);
  fftw_execute(ifft_of_rows4);
  for (int k = 0; k < M_ * M_ * M_; k++)
  {
    (*d2)[k] = conj((*d2)[k]);
    (*d4)[k] = conj((*d4)[k]);
    (*d2)[k] = (*d1)[k] + (*d2)[k];
    (*d4)[k] = (*d3)[k] + (*d4)[k];
  }
  // start here we need to transpose the arrays first
  // We need cannot do a pure 1D decomposition of 3D FFT
  // transpose
  for (int j = 0; j < M_; j++)
  {
    for (int k = 0; k < M_; k++)
    {
      for (int l = 0; l < M_; l++)
      {
        ((*d1)[l + (M_) * (k + (M_)*j)]).real(((*d2)[k + (M_) * (l + (M_)*j)]).real());
        ((*d1)[l + (M_) * (k + (M_)*j)]).imag(((*d2)[k + (M_) * (l + (M_)*j)]).imag());
        ((*d3)[l + (M_) * (k + (M_)*j)]).real(((*d4)[k + (M_) * (l + (M_)*j)]).real());
        ((*d3)[l + (M_) * (k + (M_)*j)]).imag(((*d4)[k + (M_) * (l + (M_)*j)]).imag());
      }
    }
  }
  fftw_execute(ifft_of_collums1);
  fftw_execute(ifft_of_collums2);
  for (int j = 0; j < M_; j++)
  {
    for (int k = 0; k < M_; k++)
    {
      for (int l = 0; l < M_; l++)
      {
        ((*d1)[l + (M_) * (k + (M_)*j)]).real(((*d5)[k + (M_) * (l + (M_)*j)]).real());
        ((*d1)[l + (M_) * (k + (M_)*j)]).imag(((*d5)[k + (M_) * (l + (M_)*j)]).imag());
        ((*d3)[l + (M_) * (k + (M_)*j)]).real(((*d6)[k + (M_) * (l + (M_)*j)]).real());
        ((*d3)[l + (M_) * (k + (M_)*j)]).imag(((*d6)[k + (M_) * (l + (M_)*j)]).imag());
      }
    }
  }
  // and back here

  for (int k = 0; k < M_ * M_ * M_; k++)
  {
    (*d6)[k] = conj((*d3)[k]);
    (*d5)[k] = (*d1)[k] + (*d6)[k];
  }
  fftw_execute(ifft_of_rank);
  for (int i = 0; i < M_ * M_ * M_; i++)
  {
    values_[i] = ((*d7)[i]).real();
  }
  fftw_destroy_plan(ifft_of_rows1);
  fftw_destroy_plan(ifft_of_rows2);
  fftw_destroy_plan(ifft_of_rows3);
  fftw_destroy_plan(ifft_of_rows4);
  fftw_destroy_plan(ifft_of_collums1);
  fftw_destroy_plan(ifft_of_collums2);
  fftw_destroy_plan(ifft_of_rank);
}
double UQ::RandomFieldSpectral::SimGaussRandomFieldCOS3D(double x, double y, double z)
{
  double result = 0;
  for (int j = 0; j < N_; j++)
  {
    for (int k = 0; k < N_; k++)
    {
      for (int l = 0; l < N_; l++)
      {
        // set first elements to zero
        // if(k==0||j==0||l==0||j>(N_-2)||k>(N_-2)||l>(N_-2))
        result += sqrt(2 * discrete_PSD_[l + N_ * (k + N_ * j)] * pow(dkappa_, 3)) *
                  (cos((l)*dkappa_ * x + (k)*dkappa_ * y + (j)*dkappa_ * z +
                       Phi_0_[l + N_ * (k + N_ * j)]) +
                      cos((l)*dkappa_ * x + (k)*dkappa_ * y - (j)*dkappa_ * z +
                          Phi_1_[l + N_ * (k + N_ * j)]) +
                      cos((l)*dkappa_ * x - (k)*dkappa_ * y - (j)*dkappa_ * z +
                          Phi_2_[l + N_ * (k + N_ * j)]) +
                      cos((l)*dkappa_ * x - (k)*dkappa_ * y + (j)*dkappa_ * z +
                          Phi_3_[l + N_ * (k + N_ * j)]));
      }
    }
  }
  return sqrt(2) * result;
}


/*double UQ::RandomFieldSpectral::EvalFieldAtLocation(std::vector<double> location, bool
writetofile, bool output)
{
  // manage the two different variants for evalutation in here so that it cannot be seen from the
outside
  // and so that we can call the same function with the same syntax
 if (UseFFT_)
  {
    int index_x;
    int index_y;
    int index_z;

    // Compute indices
    index_x=int(floor((location[0]-bb_min_[0])/dx_));
    index_y=int(floor((location[1]-bb_min_[1])/dx_));
    // HACK for 2D art_aorta_case SET z to y
     if (myrank_ == 0&& output && dim_==2 )
     {
       IO::cout<< "hack in use" << IO::endl;
     }
     if (dim_==2)
       index_y=int(floor((location[2]-bb_min_[2])/dx_));
       index_z=int(floor((location[2]-bb_min_[2])/dx_));
    // check index
    if (index_x>M_||index_y>M_||index_z>M_)
      dserror("Index out of bounds");
    if (writetofile && myrank_==0 )
    {
      std::ofstream File;
      File.open("RFatPoint.txt",std::ios::app);
      // use at() to get an error massage just in case
      File << std::setprecision (9) << values_[index_x+M_*index_y]<< std::endl;
      File.close();
    }
    if (dim_==2)
    {
      return values_[index_x+M_*index_y];
    }
    else
      return values_[index_x+M_*(index_y+M_*index_z)];
  }
 else
 {
   if(dim_==3)
   {
     double value;
     value=SimGaussRandomFieldCOS3D(location[0], location[1], location[2]);
     TranslateToNonGaussian( &value);
     if (writetofile && myrank_==0 )
         {
           std::ofstream File;
           File.open("RFatPoint.txt",std::ios::app);
           File << std::setprecision (9) << value<< std::endl;
           File.close();
         }
     return value;
   }
   else
   {
     dserror("Computation using Cos series only for dim = 3" );
     return -1;
   }

 }

}*/

double UQ::RandomFieldSpectral::EvalFieldAtLocation(const std::vector<double> location,
    const double paracont_parameter, const bool writetofile, const bool output)
{
  // manage the two different variants for evalutation in here so that it cannot be seen from the
  // outside and so that we can call the same function with the same syntax

  double temp_rf_val = -1.0;

  if (UseFFT_)
  {
    int index_x;
    int index_y;
    int index_z;

    // Compute indices
    index_x = int(floor((location[0] - bb_min_[0]) / dx_));
    index_y = int(floor((location[1] - bb_min_[1]) / dx_));
    // special solution for 2D art_aorta_case SET z to y
    if (myrank_ == 0 && output && dim_ == 2)
    {
      IO::cout << "hack in use" << IO::endl;
    }
    if (dim_ == 2) index_y = int(floor((location[2] - bb_min_[2]) / dx_));
    index_z = int(floor((location[2] - bb_min_[2]) / dx_));
    // check index
    if (index_x > M_ || index_y > M_ || index_z > M_) dserror("Index out of bounds");
    if (writetofile && myrank_ == 0)
    {
      std::ofstream File;
      File.open("RFatPoint.txt", std::ios::app);
      // use at() to get an error massage just in case
      File << std::setprecision(9) << values_[index_x + M_ * index_y] << std::endl;
      File.close();
    }
    if (dim_ == 2)
    {
      temp_rf_val = det_value_paracont_ * (1.0 - paracont_parameter) +
                    (paracont_parameter)*values_[index_x + M_ * index_y];
      // return values_[index_x+M_*index_y];
    }
    else
    {
      // return values_[index_x+M_*(index_y+M_*index_z)];
      temp_rf_val = det_value_paracont_ * (1.0 - paracont_parameter) +
                    (paracont_parameter)*values_[index_x + M_ * (index_y + M_ * index_z)];
    }
  }
  else
  {
    if (dim_ == 3)
    {
      double value;
      value = SimGaussRandomFieldCOS3D(location[0], location[1], location[2]);
      TranslateToNonGaussian(&value);
      if (writetofile && myrank_ == 0)
      {
        std::ofstream File;
        File.open("RFatPoint.txt", std::ios::app);
        File << std::setprecision(9) << value << std::endl;
        File.close();
      }
      // return value;
      temp_rf_val = det_value_paracont_ * (1.0 - paracont_parameter) + (paracont_parameter)*value;
    }
    else
    {
      dserror("Computation using Cos series only for dim = 3");
    }
  }
  // quick check if the value if outside of specified bounds
  if (is_bounded_)
  {
    if (temp_rf_val > rf_upper_bound_) temp_rf_val = rf_upper_bound_;
    if (temp_rf_val < rf_lower_bound_) temp_rf_val = rf_lower_bound_;
  }
  return temp_rf_val;
}


// Translate Gaussian to nonGaussian process based on Mircea Grigoriu's translation process
// theory
void UQ::RandomFieldSpectral::TranslateToNonGaussian()
{
  double dim = dim_;
  // check wether pdf is gaussian
  switch (marginal_pdf_)
  {
    case normal:
      // IO::cout<< RED_LIGHT << "WARNING: Target marginal PDF is gaussian so nothing to do here"<<
      // END_COLOR << std::endl;
      break;

    case beta:
    {
      dserror("Beta distribution not supported yet");
      //      for(int i=0;i<(pow(M_,dim));i++)
      //      {
      //        values_[i]=boost::math::quantile(my_beta,boost::math::cdf(my_norm,
      //        values_[i]))*distribution_params_[3]+distribution_params_[2];;
      //      }
    }
    break;

    case lognormal:
    {
      // This works if we use FFT
      boost::accumulators::accumulator_set<double, stats<boost::accumulators::tag::variance>> acc;
      for (int i = 0; i < (pow(M_, dim)); i++)
      {
        acc(values_[i]);
      }
      // normal_distribution<>  my_norm2(0,sqrt(variance(acc)));
      // End of FFT

      // normal_distribution<>  my_norm2(0,sigma_ul_g_cur_it_);
      boost::math::normal_distribution<double> my_norm2(0, sigma_ul_g_cur_it_);
      boost::math::lognormal_distribution<double> my_lognorm(
          distribution_params_[0], distribution_params_[1]);


      IO::cout << "translating with sigma_ul_g_cur_it_= " << sigma_ul_g_cur_it_ << IO::endl;
      IO::cout << "distribution_params_[0]" << distribution_params_[0] << IO::endl;
      IO::cout << "distribution_params_[1]" << distribution_params_[1] << IO::endl;

      // actual translation see Grigoriu1995 for details
      for (int i = 0; i < (pow(M_, dim)); i++)
      {
        values_[i] = boost::math::quantile(my_lognorm, boost::math::cdf(my_norm2, values_[i]));
      }
    }
    break;
    default:
      dserror("Only lognormal and beta distribution supported fix your input file");
      break;
  }
}

double UQ::RandomFieldSpectral::EvalFieldAtLocation(std::vector<double> location,
    double local_median, double paracont_parameter, bool writetofile, bool output)
{
  dserror(
      "EvalFieldAtLocation with local median value as argument not implemented for spectral "
      "field!");
  return -1.0;
}
// Overloaded function to translate single point only
void UQ::RandomFieldSpectral::TranslateToNonGaussian(double* value)
{
  // check wether pdf is gaussian
  switch (marginal_pdf_)
  {
    case normal:
      // IO::cout<< RED_LIGHT << "WARNING: Target marginal PDF is gaussian so nothing to do here"<<
      // END_COLOR << IO::endl;
      break;
    case beta:
    {
      dserror("fix this function");
    }
    break;

    case lognormal:
    {
      // estimate variance from PSD
      boost::math::normal_distribution<double> my_norm2(0, sigma_ul_g_cur_it_);
      boost::math::lognormal_distribution<double> my_lognorm(
          distribution_params_[0], distribution_params_[1]);
      // IO::cout<< "distribution_params_[0]" << distribution_params_[0] << IO::endl;
      // IO::cout<< "distribution_params_[1]" << distribution_params_[1] << IO::endl;
      *value = boost::math::quantile(my_lognorm, boost::math::cdf(my_norm2, *value));
    }
    break;
    default:
      dserror("Only lognormal and beta distribution supported fix your input file");
      break;
  }
}
// Transform PSD of underlying gauusian process
void UQ::RandomFieldSpectral::SpectralMatching()
{
  // error to target psd init with 100 %
  double psd_error = 100;
  double error_numerator;
  double error_denominator;
  // Target PSD of non gassian field
  std::vector<double> PSD_ng_target;
  std::vector<double> PSD_ng(N_ * N_, 0.0);
  std::vector<double> PSD_ul_g(N_ * N_, 0.0);
  std::vector<double> rho(2 * N_ * 2 * N_, 0.0);

  PSD_ng_target = discrete_PSD_;

  PSD_ng_target[0] = 0.0;


  Teuchos::RCP<Teuchos::Array<std::complex<double>>> autocorr =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> almost_autocorr =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> autocorr_ng =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2));

  Teuchos::RCP<Teuchos::Array<std::complex<double>>> PSD_ul_g_complex =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> almost_PSD_ng_complex =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> PSD_ng_complex =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2));

  for (int j = 0; j < N_ * 2; j++)
  {
    for (int k = 0; k < N_ * 2; k++)
    {
      // sort entries ro w major style
      if (j > (N_ - 1) || k > (N_ - 1))
      {
        ((*PSD_ul_g_complex)[k + N_ * 2 * j]).real(0.0);
        ((*PSD_ul_g_complex)[k + N_ * 2 * j]).imag(0.0);
      }
      else
      {
        ((*PSD_ul_g_complex)[k + N_ * 2 * j]).real(discrete_PSD_[k + j * N_]);
        ((*PSD_ul_g_complex)[k + N_ * 2 * j]).imag(0.0);
        if (k == 0 || j == 0)
        {
          // we have 0.5 already in computation of discrete psd
          // r//eal((*PSD_ul_g_complex)[k+N_*2*j])=std::real((*PSD_ul_g_complex)[k+N_*2*j])*0.5;
          // std::imag((*PSD_ul_g_complex)[k+N_*2*j])=0.0;
        }
        else if (k == 0 && j == 0)
        {
          // std::real((*PSD_ul_g_complex)[k+N_*2*j])=0;
          // std::real((*PSD_ul_g_complex)[k+N_*2*j])=std::real((*PSD_ul_g_complex)[k+N_*2*j])*0.25;
          // std::imag((*PSD_ul_g_complex)[k+N_*2*j])=0;
        }
      }
    }
  }
  int rank =
      1; /* not 2: we are computing 1d transforms int n[] = {1024}; 1d transforms of length 10 */
  int N_fftw = 2 * N_;
  int howmany = 2 * N_;  // same here
  int idist = 2 * N_;
  int odist = 2 * N_;
  int istride = 1;
  int ostride = 1; /* distance between two elements in the same row/column/rank */


  fftw_plan ifft_of_rows_of_psd;
  fftw_plan ifft_of_columns_of_psd;
  // ifft for autocorr
  fftw_plan ifft_of_rows_of_autocorr_ng;
  fftw_plan ifft_of_columns_of_autocorr_ng;

  ifft_of_rows_of_psd = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*PSD_ul_g_complex)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*almost_autocorr)[0]))), NULL, ostride, odist,
      FFTW_BACKWARD, FFTW_ESTIMATE);

  ifft_of_rows_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*almost_PSD_ng_complex)[0]))), NULL, ostride, odist,
      FFTW_FORWARD, FFTW_ESTIMATE);


  istride = 2 * N_;
  ostride = 2 * N_;
  idist = 1;
  odist = 1;
  ifft_of_columns_of_psd = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*almost_autocorr)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))), NULL, ostride, odist, FFTW_BACKWARD,
      FFTW_ESTIMATE);
  ifft_of_columns_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*almost_PSD_ng_complex)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))), NULL, ostride, odist,
      FFTW_FORWARD, FFTW_ESTIMATE);

  sigma_ul_g_cur_it_ = sigma_0_;


  // while error > 0.5%
  int i = 0;
  do
  {
    if (i != 0)  // set new psd_ul_g
    {
      for (int j = 0; j < N_ * 2; j++)
      {
        for (int k = 0; k < N_ * 2; k++)
        {
          if (j > (N_ - 1) || k > (N_ - 1))
          {
            ((*PSD_ul_g_complex)[k + N_ * 2 * j]).real(0.0);
            ((*PSD_ul_g_complex)[k + N_ * 2 * j]).imag(0.0);
          }
          else
          {
            ((*PSD_ul_g_complex)[k + N_ * 2 * j]).real(PSD_ul_g[k + j * N_]);
            ((*PSD_ul_g_complex)[k + N_ * 2 * j]).imag(0.0);
            if (k == 0 || j == 0)
            {
              // was working when commented out
              // std::real((*PSD_ul_g_complex)[k+N_*2*j])=std::real((*PSD_ul_g_complex)[k+N_*2*j])*0.5;
              // factor 0.5 taken care of by updating formular
            }
          }
        }
      }
      double sigma_ul_g_cur_it_helper_ = 0;
      for (int g = 0; g < 4 * N_ * N_; g++)
      {
        sigma_ul_g_cur_it_helper_ += ((*PSD_ul_g_complex)[g]).real() * pow(dkappa_, 2);
      }
      sigma_ul_g_cur_it_helper_ =
          sigma_ul_g_cur_it_helper_ - ((*PSD_ul_g_complex)[0]).real() * pow(dkappa_, 2) * 1.0;
      sigma_ul_g_cur_it_ = sqrt(4 * sigma_ul_g_cur_it_helper_);
      IO::cout << "Sigma of PSD_UL_G" << sigma_ul_g_cur_it_ << IO::endl;
    }
    fftw_execute(ifft_of_rows_of_psd);
    fftw_execute(ifft_of_columns_of_psd);
    double scaling_fac = dkappa_ * N_;

    // loop over vectorlength
    for (int k = 0; k < 2 * N_ * 2 * N_; k++)
    {
      rho[k] = ((*autocorr)[k]).real() * 2 * pow((scaling_fac), 2) /
               (2 * N_ * 2 * N_ * (pow(sigma_0_, 2)));
      // lets go for +_ 3 * sigma here
      ((*autocorr_ng)[k])
          .real(Integrate(-3 * sigma_ul_g_cur_it_, 3 * sigma_ul_g_cur_it_, -3 * sigma_ul_g_cur_it_,
              3 * sigma_ul_g_cur_it_, rho[k]));
    }
    fftw_execute(ifft_of_rows_of_autocorr_ng);
    fftw_execute(ifft_of_columns_of_autocorr_ng);

    for (int j = 0; j < N_ * 2; j++)
    {
      for (int k = 0; k < N_ * 2; k++)
      {
        // sort entries ro w major style
        // set first elements to zero
        // if(k==0||j==0||j>(N_-1)||k>(N_-1))
        if (j > (N_ - 1) || k > (N_ - 1))
        {
        }
        // else if (j==0&& k== 0)
        // PSD_ng[k+j*N_]=0.0;
        else
        {
          PSD_ng[k + j * N_] = ((*PSD_ng_complex)[k + j * 2 * N_]).real() / (pow(scaling_fac, 2));
          PSD_ul_g[k + j * N_] = ((*PSD_ul_g_complex)[k + j * 2 * N_]).real();
        }
      }
    }

    PSD_ng[0] = 0.0;

    for (int k = 0; k < N_ * N_; k++)
    {
      if (PSD_ng[k] > 10e-10) PSD_ul_g[k] = pow(PSD_ng_target[k] / PSD_ng[k], 1.4) * PSD_ul_g[k];
      // do not set to zero because if once zero you'll never get it non-zero again
      else
        PSD_ul_g[k] = 10e-10;
    }
    // compute error based on equation(19) from shield2011
    error_numerator = 0.0;
    error_denominator = 0.0;
    for (int g = 0; g < N_ * N_; g++)
    {
      error_numerator += pow((PSD_ng[g] - PSD_ng_target[g]), 2);
      error_denominator += pow((PSD_ng_target[g]), 2);
    }
    psd_error = 100 * sqrt(error_numerator / error_denominator);
    if (myrank_ == 0) IO::cout << "Error to target PSD: " << psd_error << IO::endl;
    // increase counter
    i++;
  }
  // set error threshold for spectral matching to 0.5 %
  while (psd_error > 0.1);



  // free memory
  fftw_destroy_plan(ifft_of_columns_of_psd);
  fftw_destroy_plan(ifft_of_rows_of_psd);
  fftw_destroy_plan(ifft_of_rows_of_autocorr_ng);
  fftw_destroy_plan(ifft_of_columns_of_autocorr_ng);


  // Write PSD_ul_g PSD_ng_target and PSD_ng to a file
  // Dimension is 128* 128

  for (int h = 0; h < N_ * N_; h++)
  {
    if (PSD_ng[h] > 10e-10)  // change that
      discrete_PSD_[h] = PSD_ul_g[h];
    else
      // remove all the very small entries to get rid of the wiggles
      discrete_PSD_[h] = 0.0;
  }
  // Write to file
  /*if (myrank_ == 0 && !reduced_output_)
    {
      std::ofstream File;
      File.open("DiscretePSDTranslated.txt",std::ios::out);
      int size = int (pow(N_,2.0));
      for(int i=0;i<size;i++)
      {
        File << discrete_PSD_[i]<< std::endl;
      }
      File.close();
      std::ofstream File2;
      File2.open("PSD_ng.txt",std::ios::out);
        for(int i=0;i<size;i++)
        {
          File2 << PSD_ng[i]<< std::endl;
        }
        File2.close();
    }*/
  IO::cout << "Spectral Matching done " << IO::endl;
}
// Transform PSD of underlying gauusian process
void UQ::RandomFieldSpectral::SpectralMatching3D()
{
  // error to target psd init with 100 %
  double psd_error = 100;
  double error_numerator;
  double error_denominator;
  // Target PSD of non gassian field
  std::vector<double> PSD_ng_target;
  std::vector<double> PSD_ng(N_ * N_ * N_, 0.0);
  std::vector<double> PSD_ul_g(N_ * N_ * N_, 0.0);
  std::vector<double> rho(2 * N_ * 2 * N_ * 2 * N_, 0.0);

  PSD_ng_target = discrete_PSD_;

  PSD_ng_target[0] = 0.0;


  Teuchos::RCP<Teuchos::Array<std::complex<double>>> autocorr =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> almost_autocorr =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> almost_autocorr2 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> autocorr_ng =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));

  Teuchos::RCP<Teuchos::Array<std::complex<double>>> PSD_ul_g_complex =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> almost_PSD_ng_complex =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> almost_PSD_ng_complex2 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> PSD_ng_complex =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> temp =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));

  for (int j = 0; j < N_ * 2; j++)
  {
    for (int k = 0; k < N_ * 2; k++)
    {
      for (int l = 0; l < N_ * 2; l++)
      {
        // sort entries ro w major style
        if (j > (N_ - 1) || k > (N_ - 1) || l > (N_ - 1))
        {
          ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).real(0.0);
          ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).imag(0.0);
        }
        else
        {
          ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)])
              .real(discrete_PSD_[l + (N_) * (k + (N_)*j)]);
          ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).imag(0.0);
        }
      }
    }
  }

  int rank =
      1; /* not 3: we are computing 1d transforms int n[] = {1024} 1d transforms of length 10 */
  int N_fftw = 2 * N_;
  int howmany = 2 * N_ * 2 * N_;  // same here
  int idist = 2 * N_;  // the distance in memory between the first element  of the first array and
                       // the first element of the second array */
  int odist = 2 * N_;
  int istride = 1;
  int ostride = 1; /* distance between two elements in the same row/column/rank */

  fftw_plan ifft_of_rows_of_psd;
  fftw_plan ifft_of_columns_of_psd;
  fftw_plan ifft_of_rank_of_psd;
  // ifft for autocorr
  fftw_plan ifft_of_rows_of_autocorr_ng;
  fftw_plan ifft_of_columns_of_autocorr_ng;
  fftw_plan ifft_of_rank_of_autocorr_ng;


  ifft_of_rows_of_psd = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*PSD_ul_g_complex)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))), NULL, ostride, odist, FFTW_BACKWARD,
      FFTW_ESTIMATE);


  ifft_of_rows_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))), NULL, ostride, odist,
      FFTW_BACKWARD, FFTW_ESTIMATE);

  // now the same for collumns
  idist = 2 * N_;
  odist = 2 * N_;  // the distance in memory between the first element  of the first array and the
                   // first element of the second array */
  istride = 1;     /* distance between two elements in the same row/column/rank */
  ostride = 1;     /* distance between two elements in the same row/column/rank */
  howmany = 2 * N_ * 2 * N_;  // just 2*N_ we call this multiple times
  rank = 1; /* not 3: we are computing 1d transforms int n[] = {1024} 1d transforms of length 10 */


  ifft_of_columns_of_psd =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*temp)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*almost_autocorr2)[0]))), NULL,
          ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);

  ifft_of_columns_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*almost_PSD_ng_complex2)[0]))), NULL, ostride, odist,
      FFTW_BACKWARD, FFTW_ESTIMATE);

  // and now the rank (aka third dim of the array)
  howmany = 2 * N_ * 2 * N_;
  idist = 1;  // the distance in memory between the first element  of the first array and the first
              // element of the second array */
  odist = 1;
  istride = 2 * N_ * 2 * N_;
  ostride = 2 * N_ * 2 * N_; /* distance between two elements in the same row/column/rank */
  ifft_of_rank_of_psd = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*almost_autocorr)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))), NULL, ostride, odist, FFTW_BACKWARD,
      FFTW_ESTIMATE);

  ifft_of_rank_of_autocorr_ng = fftw_plan_many_dft(rank, &N_fftw, howmany,
      (reinterpret_cast<fftw_complex*>(&((*almost_PSD_ng_complex)[0]))), NULL, istride, idist,
      (reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))), NULL, ostride, odist,
      FFTW_BACKWARD, FFTW_ESTIMATE);

  sigma_ul_g_cur_it_ = sigma_0_;
  // while error > 0.5%
  int i = 0;
  do
  {
    if (i != 0)  // set new psd_ul_g
    {
      for (int j = 0; j < N_ * 2; j++)
      {
        for (int k = 0; k < N_ * 2; k++)
        {
          for (int l = 0; l < N_ * 2; l++)
          {
            if (j > (N_ - 1) || k > (N_ - 1) || l > (N_ - 1))
            {
              ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).real(0.0);
              ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).imag(0.0);
            }
            else
            {
              ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)])
                  .real(PSD_ul_g[l + (N_) * (k + (N_)*j)]);
              ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).imag(0.0);
            }
          }
        }
      }
      double sigma_ul_g_cur_it_helper_ = 0;
      for (int g = 0; g < 8 * N_ * N_ * N_; g++)
      {
        sigma_ul_g_cur_it_helper_ += ((*PSD_ul_g_complex)[g]).real() * pow(dkappa_, 3);
      }
      sigma_ul_g_cur_it_ = sqrt(8 * sigma_ul_g_cur_it_helper_);
      IO::cout << "Sigma of PSD_UL_G" << sigma_ul_g_cur_it_ << IO::endl;
    }
    IO::cout << "Sigma of PSD_UL_G" << sigma_ul_g_cur_it_ << IO::endl;

    fftw_execute(ifft_of_rows_of_psd);

    for (int j = 0; j < N_ * 2; j++)
    {
      for (int k = 0; k < N_ * 2; k++)
      {
        for (int l = 0; l < N_ * 2; l++)
        {
          ((*temp)[l + (2 * N_) * (k + (2 * N_) * j)])
              .real(((*autocorr)[k + (2 * N_) * (l + (2 * N_) * j)]).real());
          ((*temp)[l + (2 * N_) * (k + (2 * N_) * j)])
              .imag(((*autocorr)[k + (2 * N_) * (l + (2 * N_) * j)]).imag());
        }
      }
    }
    fftw_execute(ifft_of_columns_of_psd);
    // and transpose back
    // transpose
    for (int j = 0; j < N_ * 2; j++)
    {
      for (int k = 0; k < N_ * 2; k++)
      {
        for (int l = 0; l < N_ * 2; l++)
        {
          ((*almost_autocorr)[l + (2 * N_) * (k + (2 * N_) * j)])
              .real(((*almost_autocorr2)[k + (2 * N_) * (l + (2 * N_) * j)]).real());
          // if you forget the std::imaginary part it will mess with your scaling
          ((*almost_autocorr)[l + (2 * N_) * (k + (2 * N_) * j)])
              .imag(((*almost_autocorr2)[k + (2 * N_) * (l + (2 * N_) * j)]).imag());
        }
      }
    }


    fftw_execute(ifft_of_rank_of_psd);
    double scaling_fac = dkappa_ * N_;

    for (int k = 0; k < 2 * N_ * 2 * N_ * 2 * N_; k++)
    {
      // scaling with sigma to get autocorrelation function
      // Factor 2 to here at the back is essential (altough not quite sure were it comes from)
      rho[k] = ((*autocorr)[k]).real() * 2 * pow((scaling_fac), 3) /
               (2 * N_ * 2 * N_ * 2 * N_ * (pow(sigma_0_, 2)));

      ((*autocorr_ng)[k])
          .real(Integrate(-3 * sigma_ul_g_cur_it_, 3 * sigma_ul_g_cur_it_, -3 * sigma_ul_g_cur_it_,
              3 * sigma_ul_g_cur_it_, rho[k]));
      ((*autocorr_ng)[k]).imag(0.0);
      // The followign lines are good for testing if the FFT works correctly
      // rho[k] = std::real((*autocorr)[k])*pow((scaling_fac),3)/(2*N_*2*N_*2*N_);
    }
    fftw_execute(ifft_of_rows_of_autocorr_ng);
    // transpose
    for (int j = 0; j < N_ * 2; j++)
    {
      for (int k = 0; k < N_ * 2; k++)
      {
        for (int l = 0; l < N_ * 2; l++)
        {
          ((*autocorr_ng)[l + (2 * N_) * (k + (2 * N_) * j)])
              .real(((*PSD_ng_complex)[k + (2 * N_) * (l + (2 * N_) * j)]).real());
          ((*autocorr_ng)[l + (2 * N_) * (k + (2 * N_) * j)])
              .imag(((*PSD_ng_complex)[k + (2 * N_) * (l + (2 * N_) * j)]).imag());
        }
      }
    }
    fftw_execute(ifft_of_columns_of_autocorr_ng);
    // transpose back
    for (int j = 0; j < N_ * 2; j++)
    {
      for (int k = 0; k < N_ * 2; k++)
      {
        for (int l = 0; l < N_ * 2; l++)
        {
          ((*almost_PSD_ng_complex)[l + (2 * N_) * (k + (2 * N_) * j)])
              .real(((*almost_PSD_ng_complex2)[k + (2 * N_) * (l + (2 * N_) * j)]).real());
          ((*almost_PSD_ng_complex)[l + (2 * N_) * (k + (2 * N_) * j)])
              .imag(((*almost_PSD_ng_complex2)[k + (2 * N_) * (l + (2 * N_) * j)]).imag());
        }
      }
    }


    fftw_execute(ifft_of_rank_of_autocorr_ng);
    // IO::cout<< "PSD after fft ( 10 10 10) "<<
    // std::real((*PSD_ng_complex)[10+(2*N_)*(10+(2*N_)*2)])/pow((scaling_fac),3)*2<< "std::imag "<<
    // std::imag((*PSD_ng_complex)[10+(2*N_)*(10+(2*N_)*2)])/pow((scaling_fac),3)*8.  << IO::endl;
    // IO::cout<< "PSD after fft ( 16644) "<<
    // std::real((*PSD_ng_complex)[16644])/pow((scaling_fac),3)*2 << IO::endl;
    for (int j = 0; j < N_ * 2; j++)
    {
      for (int k = 0; k < N_ * 2; k++)
      {
        for (int l = 0; l < N_ * 2; l++)
        {
          if (j > (N_ - 1) || k > (N_ - 1) || l > (N_ - 1))
          {
          }
          else
          {
            PSD_ng[l + (N_) * (k + (N_)*j)] =
                ((*PSD_ng_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).real() /
                (pow(scaling_fac, 3));
            PSD_ul_g[l + (N_) * (k + (N_)*j)] =
                ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).real();
          }
        }
      }
    }

    PSD_ng[0] = 0.0;
    for (int k = 0; k < N_ * N_ * N_; k++)
    {
      if (PSD_ng[k] > 10e-10)
      {
        PSD_ul_g[k] = pow(PSD_ng_target[k] / PSD_ng[k], 1.4) * PSD_ul_g[k];
      }
      // do not set to zero because if once zero you'll never get it non-zero again
      else
      {
        PSD_ul_g[k] = 10e-10;
      }
    }

    // compute error based on equation(19) from shield2011
    error_numerator = 0.0;
    error_denominator = 0.0;
    for (int g = 0; g < N_ * N_ * N_; g++)
    {
      error_numerator += pow((PSD_ng[g] - PSD_ng_target[g]), 2);
      error_denominator += pow((PSD_ng_target[g]), 2);
    }
    psd_error = 100 * sqrt(error_numerator / error_denominator);
    if (myrank_ == 0) IO::cout << "Error to target PSD: " << psd_error << IO::endl;
    i++;
  }
  // set error threshold for spectral matching to 0.5 %
  while (psd_error > 0.5);
  // free memory
  fftw_destroy_plan(ifft_of_rows_of_psd);
  fftw_destroy_plan(ifft_of_columns_of_psd);
  fftw_destroy_plan(ifft_of_rank_of_psd);
  fftw_destroy_plan(ifft_of_rows_of_autocorr_ng);
  fftw_destroy_plan(ifft_of_columns_of_autocorr_ng);
  fftw_destroy_plan(ifft_of_rank_of_autocorr_ng);

  for (int h = 0; h < N_ * N_ * N_; h++)
  {
    if (PSD_ng[h] > 10e-10)  // change that
      discrete_PSD_[h] = PSD_ul_g[h];
    else
      // remove all the very small entries to get rid of the wiggles
      discrete_PSD_[h] = 0.0;
  }
}

// Routine to calculate
// Transform PSD of underlying gauusian process using direct 3DFFT routine provided by FFTW
void UQ::RandomFieldSpectral::SpectralMatching3D3D()
{
  // error to target psd init with 100 %
  double psd_error = 100;
  double error_numerator;
  double error_denominator;
  // Target PSD of non gassian field
  std::vector<double> PSD_ng_target;
  std::vector<double> PSD_ng(N_ * N_ * N_, 0.0);
  std::vector<double> PSD_ul_g(N_ * N_ * N_, 0.0);
  std::vector<double> rho(2 * N_ * 2 * N_ * 2 * N_, 0.0);

  PSD_ng_target = discrete_PSD_;
  // do i need to set this zero??
  PSD_ng_target[0] = 0.0;
  // calc sigma form discrete PSD

  Teuchos::RCP<Teuchos::Array<std::complex<double>>> autocorr =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> autocorr_ng =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));

  Teuchos::RCP<Teuchos::Array<std::complex<double>>> PSD_ul_g_complex =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> PSD_ng_complex =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(N_ * 2 * N_ * 2 * N_ * 2, 0.0));

  for (int j = 0; j < N_ * 2; j++)
  {
    for (int k = 0; k < N_ * 2; k++)
    {
      for (int l = 0; l < N_ * 2; l++)
      {
        // sort entries ro w major style
        if (j > (N_ - 1) || k > (N_ - 1) || l > (N_ - 1))
        {
          // l+M_*(k+M_*j)
          ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).real(0.0);
          ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).imag(0.0);
        }
        else
        {
          ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)])
              .real(discrete_PSD_[l + (N_) * (k + (N_)*j)]);
          ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).imag(0.0);
        }
      }
    }
  }
  fftw_plan fft_of_psd;
  fftw_plan fft_of_autocorr_ng;
  // sign, can be either FFTW_FORWARD (-1) or FFTW_BACKWARD (+1),
  fft_of_psd = fftw_plan_dft_3d(2 * N_, 2 * N_, 2 * N_,
      (reinterpret_cast<fftw_complex*>(&((*PSD_ul_g_complex)[0]))),
      (reinterpret_cast<fftw_complex*>(&((*autocorr)[0]))), 1, FFTW_ESTIMATE);

  fft_of_autocorr_ng = fftw_plan_dft_3d(2 * N_, 2 * N_, 2 * N_,
      (reinterpret_cast<fftw_complex*>(&((*autocorr_ng)[0]))),
      (reinterpret_cast<fftw_complex*>(&((*PSD_ng_complex)[0]))), 1, FFTW_ESTIMATE);
  sigma_ul_g_cur_it_ = sigma_0_;
  // while error > 0.5%
  int i = 0;
  do
  {
    if (i != 0)  // set new psd_ul_g
    {
      for (int j = 0; j < N_ * 2; j++)
      {
        for (int k = 0; k < N_ * 2; k++)
        {
          for (int l = 0; l < N_ * 2; l++)
          {
            if (j > (N_ - 1) || k > (N_ - 1) || l > (N_ - 1))
            {
              ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).real(0.0);
              ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).imag(0.0);
            }
            else
            {
              ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)])
                  .real(PSD_ul_g[l + (N_) * (k + (N_)*j)]);
              ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).imag(0.0);
            }
          }
        }
      }
      double sigma_ul_g_cur_it_helper_ = 0;
      for (int g = 0; g < 8 * N_ * N_ * N_; g++)
      {
        sigma_ul_g_cur_it_helper_ += ((*PSD_ul_g_complex)[g]).real() * pow(dkappa_, 3);
      }
      sigma_ul_g_cur_it_helper_ =
          sigma_ul_g_cur_it_helper_ - ((*PSD_ul_g_complex)[0]).real() * pow(dkappa_, 2) * 0.0;
      sigma_ul_g_cur_it_ = sqrt(8 * sigma_ul_g_cur_it_helper_);
      IO::cout << "Sigma of PSD_UL_G" << sigma_ul_g_cur_it_ << IO::endl;
    }

    fftw_execute(fft_of_psd);
    double scaling_fac = dkappa_ * N_;

    for (int k = 0; k < 2 * N_ * 2 * N_ * 2 * N_; k++)
    {
      // Factor 2 to here at the back is essential (altough not quite sure were it comes from)
      rho[k] = ((*autocorr)[k]).real() * 2 * pow((scaling_fac), 3) /
               (2 * N_ * 2 * N_ * 2 * N_ * (pow(sigma_0_, 2)));
      ((*autocorr_ng)[k])
          .real(Integrate(-3 * sigma_ul_g_cur_it_, 3 * sigma_ul_g_cur_it_, -3 * sigma_ul_g_cur_it_,
              3 * sigma_ul_g_cur_it_, rho[k]));
      // The followign lines are good for testing if the FFT works correctly
      // rho[k] = std::real((*autocorr)[k])*pow((scaling_fac),3)/(2*N_*2*N_*2*N_);
      // rho[k] =
      // std::real((*autocorr)[k])*2*pow((scaling_fac),3)/(2*N_*2*N_*2*N_*(pow(sigma_0_,2)));
    }
    fftw_execute(fft_of_autocorr_ng);
    // IO::cout<< "PSD after fft ( 10 10 10) "<<
    // std::real((*PSD_ng_complex)[10+(2*N_)*(10+(2*N_)*2)])/pow((scaling_fac),3)*8. << "std::imag
    // "<< std::imag((*PSD_ng_complex)[10+(2*N_)*(10+(2*N_)*2)])/pow((scaling_fac),3)*8.  <<
    // IO::endl;
    // IO::cout<< "PSD after fft ( 16644) "<<
    // std::real((*PSD_ng_complex)[16644])/pow((scaling_fac),3)*8 << IO::endl;

    for (int j = 0; j < N_ * 2; j++)
    {
      for (int k = 0; k < N_ * 2; k++)
      {
        for (int l = 0; l < N_ * 2; l++)
        {
          if (j > (N_ - 1) || k > (N_ - 1) || l > (N_ - 1))
          {
          }
          else
          {
            PSD_ng[l + (N_) * (k + (N_)*j)] =
                ((*PSD_ng_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).real() /
                (pow(scaling_fac, 3));
            PSD_ul_g[l + (N_) * (k + (N_)*j)] =
                ((*PSD_ul_g_complex)[l + (2 * N_) * (k + (2 * N_) * j)]).real();
          }
        }
      }
    }

    PSD_ng[0] = 0.0;

    for (int k = 0; k < N_ * N_ * N_; k++)
    {
      if (PSD_ng[k] > 10e-10) PSD_ul_g[k] = pow(PSD_ng_target[k] / PSD_ng[k], 1.4) * PSD_ul_g[k];
      // do not set to zero because if once zero you'll never get it non-zero again
      else
      {
        PSD_ul_g[k] = 10e-10;
      }
    }
    // compute error based on equation(19) from shield2011
    error_numerator = 0.0;
    error_denominator = 0.0;
    for (int g = 0; g < N_ * N_ * N_; g++)
    {
      error_numerator += pow((PSD_ng[g] - PSD_ng_target[g]), 2);
      error_denominator += pow((PSD_ng_target[g]), 2);
    }
    psd_error = 100 * sqrt(error_numerator / error_denominator);
    if (myrank_ == 0) IO::cout << "Error to target PSD: " << psd_error << IO::endl;
    // increase counter
    i++;
  }
  // set error threshold for spectral matching to 0.5 %
  while (psd_error > 0.5);
  // free memory
  fftw_destroy_plan(fft_of_psd);
  fftw_destroy_plan(fft_of_autocorr_ng);
  for (int h = 0; h < N_ * N_ * N_; h++)
  {
    if (PSD_ng[h] > 10e-10)  // change that
      discrete_PSD_[h] = PSD_ul_g[h];
    else
      // remove all the very small entries to get rid of the wiggles
      discrete_PSD_[h] = 0.0;
  }
}
double UQ::RandomFieldSpectral::Integrate(
    double xmin, double xmax, double ymin, double ymax, double rho)
{
  // get trillios gausspoints with high order
  // Teuchos::RCP<DRT::UTILS::GaussPoints> gp = DRT::UTILS::GaussPointCache::Instance().Create(
  // DRT::Element::quad4, 30 );
  Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
      DRT::UTILS::GaussPointCache::Instance().Create(DRT::Element::quad4, 14);
  // needed for transformation in [-1,1];[-1,1] space
  double hx = abs(xmax - xmin);
  double hy = abs(ymax - ymin);
  double jdet = hx * hy / 4;
  double integral_value = 0.0;

  for (int i = 0; i < gp->NumPoints(); i++)
  {
    integral_value += gp->Weight(i) * jdet *
                      Testfunction(xmin + hx / 2 * (1 + gp->Point(i)[0]),
                          ymin + hy / 2 * (1 + gp->Point(i)[1]), rho);
  }
  return integral_value;
}

double UQ::RandomFieldSpectral::Testfunction(double argument_x, double argument_y, double rho)
{
  double result = 0.0;
  boost::math::normal_distribution<double> my_normal(0, sigma_ul_g_cur_it_);
  switch (marginal_pdf_)
  {
    case lognormal:
    {
      boost::math::lognormal_distribution<double> my_lognorm(
          distribution_params_[0], distribution_params_[1]);

      result = boost::math::quantile(my_lognorm, (boost::math::cdf(my_normal, (argument_x)))) *
               boost::math::quantile(my_lognorm, (boost::math::cdf(my_normal, argument_y))) *
               (1 / (2 * pi_ * pow(sigma_ul_g_cur_it_, 2) * sqrt(1 - pow(rho, 2)))) *
               exp(-(pow((argument_x), 2) + pow((argument_y), 2) -
                       2 * rho * (argument_x) * (argument_y)) /
                   (2 * pow(sigma_ul_g_cur_it_, 2) * (1 - pow(rho, 2))));
    }
    break;
    case beta:
    {
      dserror("fix this function");
      boost::math::beta_distribution<double> my_beta(
          distribution_params_[0], distribution_params_[1]);
      result = (boost::math::quantile(my_beta, (boost::math::cdf(my_normal, (argument_x)))) *
                       distribution_params_[3] +
                   distribution_params_[2]) *
               (boost::math::quantile(my_beta, (boost::math::cdf(my_normal, argument_y))) *
                       distribution_params_[3] +
                   distribution_params_[2]) *
               (1 / (2 * pi_ * pow(sigma_0_, 2) * sqrt(1 - pow(rho, 2)))) *
               exp(-(pow((argument_x), 2) + pow((argument_y), 2) -
                       2 * rho * (argument_x) * (argument_y)) /
                   (2 * pow(sigma_0_, 2) * (1 - pow(rho, 2))));
    }
    break;
    default:
      dserror(" Lognorm distribution supported so far fix your input file");
      break;
  }
  return result;
}

// Write Random Field to file
void UQ::RandomFieldSpectral::WriteRandomFieldToFile()
{
  if (myrank_ == 0)
  {
    std::ofstream File;
    File.open("RandomField.txt", std::ios::out);
    int size = int(pow(M_, double(dim_)));
    for (int i = 0; i < size; i++)
    {
      File << std::setprecision(10) << values_[i] << std::endl;
    }
    File.close();
  }
}
// Compute PSD from current sample
void UQ::RandomFieldSpectral::GetPSDFromSample(Teuchos::RCP<Teuchos::Array<double>> sample_psd)
{
  // check wether sample_psd has correct size
  if (sample_psd->length() != M_ * M_) dserror("Sizemismatch");

  Teuchos::RCP<Teuchos::Array<std::complex<double>>> b1 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_));

  // define complex i
  std::complex<double> i_comp(0, 1);

  for (int j = 0; j < M_ * M_; j++)
  {
    ((*b1)[j]).real(values_[j]);
    ((*b1)[j]).imag(0.0);
  }
  // allocate output arrays
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d1 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_));
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d2 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_));

  fftw_plan fft_of_rows;
  fftw_plan fft_of_collums;

  int rank = 1; /* not 2: we are computing 1d transforms */
  /* 1d transforms of length M_ */
  int N_fftw = M_;
  int howmany = M_;  // same here
  int idist = M_;
  int odist = M_;
  int istride = 1;
  int ostride = 1; /* distance between two elements in the same column */

  fft_of_rows =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*b1)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d1)[0]))), NULL, ostride,
          odist, FFTW_FORWARD, FFTW_ESTIMATE);

  istride = M_;
  ostride = M_;
  idist = 1;
  odist = 1;
  fft_of_collums =
      fftw_plan_many_dft(rank, &N_fftw, howmany, (reinterpret_cast<fftw_complex*>(&((*d1)[0]))),
          NULL, istride, idist, (reinterpret_cast<fftw_complex*>(&((*d2)[0]))), NULL, ostride,
          odist, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(fft_of_rows);
  fftw_execute(fft_of_collums);

  for (int i = 0; i < M_ * M_; i++)
  {
    (*sample_psd)[i] =
        1 / (pow(dkappa_, 2)) * 1 / (M_ * M_) * 1 / (M_ * M_) * pow(abs((*d2)[i]), 2);
  }

  fftw_destroy_plan(fft_of_rows);
  fftw_destroy_plan(fft_of_collums);
}
// Compute PSD from current sample
void UQ::RandomFieldSpectral::GetPSDFromSample3D(Teuchos::RCP<Teuchos::Array<double>> sample_psd)
{
  // check wether sample_psd has correct size
  if (sample_psd->length() != M_ * M_ * M_) dserror("Sizemismatch");

  Teuchos::RCP<Teuchos::Array<std::complex<double>>> b1 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  // allocate output arrays
  Teuchos::RCP<Teuchos::Array<std::complex<double>>> d1 =
      Teuchos::rcp(new Teuchos::Array<std::complex<double>>(M_ * M_ * M_));
  // define complex i
  std::complex<double> i_comp(0, 1);

  for (int j = 0; j < M_ * M_ * M_; j++)
  {
    ((*b1)[j]).real(values_[j]);
    ((*b1)[j]).imag(0.0);
  }
  fftw_plan fft_of_rf;
  // sign, can be either FFTW_FORWARD (-1) or FFTW_BACKWARD (+1),
  fft_of_rf = fftw_plan_dft_3d(M_, M_, M_, (reinterpret_cast<fftw_complex*>(&((*b1)[0]))),
      (reinterpret_cast<fftw_complex*>(&((*d1)[0]))), -1, FFTW_ESTIMATE);

  fftw_execute(fft_of_rf);
  // move values into class variable
  for (int i = 0; i < M_ * M_ * M_; i++)
  {
    (*sample_psd)[i] =
        1 / (pow(dkappa_, 3)) * 1 / (M_ * M_ * M_) * 1 / (M_ * M_ * M_) * pow(abs((*d1)[i]), 2);
  }

  fftw_destroy_plan(fft_of_rf);
}
// Write Random Field to file
void UQ::RandomFieldSpectral::WriteSamplePSDToFile(Teuchos::RCP<Teuchos::Array<double>> sample_psd)
{
  if (myrank_ == 0)
  {
    std::ofstream File;
    File.open("SamplePSD.txt", std::ios::out);
    int size = int(pow(M_, double(dim_)));
    for (int i = 0; i < size; i++)
    {
      File << (*sample_psd)[i] << std::endl;
    }
    File.close();
  }
}

void UQ::RandomFieldSpectral::WriteRandomVariablesToFile(std::string filename)
{
  dserror("WriteRandomVariablesToFile not implemented yet for spectral field");
}

#endif
