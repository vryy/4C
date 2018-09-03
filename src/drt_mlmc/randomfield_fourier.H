/*!----------------------------------------------------------------------
\file randomfield_fourier.h
Created on: May 15, 2014
\brief
Class for generating samples of gaussian and lognormal random fields with
a squared exponential covariance structure based on a fourier series expansion.
This particular type of expansion is taken from the PhD Thesis from Tamellini
  "Polynomial approximation of PDEs with stochastic coefficients"

 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/

#ifdef HAVE_FFTW
#ifndef RANDOMFIELDFOURIER_H_
#define RANDOMFIELDFOURIER_H_

#include "randomfield.H"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_RCP.hpp>

// include random number libraries
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>

// random generator type
typedef boost::mt19937 randnumgen;


namespace UQ
{
  class RandomFieldFourier : public RandomField
  {
   public:
    //! \brief Standard constructor
    RandomFieldFourier(unsigned int seed, Teuchos::RCP<DRT::Discretization> discret,
        const Teuchos::ParameterList& rfp);

    //! \brief generates new realization of random field
    void CreateNewSample(unsigned int seed);

    void SayHello();

    //! \brief get value of random field at a specific location
    double EvalFieldAtLocation(const std::vector<double> location, const double paracont_parameter,
        const bool writetofile, const bool output);

    //! \brief get value of random field at a specific location
    double EvalFieldAtLocation(const std::vector<double> location, const double paracont_parameter,
        const double local_mean, const bool writetofile, const bool output);

    //! \brief return dimension of random field (either 2 or 3)
    int Dimension() { return dim_; }

   private:
    void CreateNewPhaseAngles(unsigned int seed);

    //! brief compute Eigenvalues
    double ComputeLambda3D(int k);

    //! \brief write random field to file
    void WriteRandomVariablesToFile(std::string filename);

   protected:
   private:
    //! \brief dimension of stochastic field
    int dim_;

    //! properties of first order probability distribution
    double mean_;
    double sigma_0_;

    //! brief size/length of random field in each dim (periodicity)
    double L_;

    //! \brief correlation length of random field descibing the decay of the correlation function
    // note this should be smaller than 0.35*L_
    double d_;

    //! brief number of terms in fourier series per dimension
    int m_;

    //! brief number of terms in fourier series after truncation (worst case m_^dim_)
    int num_k_;

    //! \brief vector to store the relevant terms of the fourier series expansion
    Teuchos::RCP<std::vector<std::vector<int>>> kb_;
    // Teuchos::RCP<Epetra_MultiVector> kb_;

    // brief vector to store the
    std::vector<double> omega_k_;

    //! \brief a random number generator object (boost)
    randnumgen randomnumbergen_;

    //! \brief generates Gaussian random numbers with zero mean and unit standard deviation
    Teuchos::RCP<boost::variate_generator<randnumgen&, boost::normal_distribution<>>> normalgen_;

    //! \brief vectors to store the random amplitudes
    Teuchos::RCP<std::vector<std::vector<double>>> xi_;

    //! \brief/ seed for random number generator
    unsigned int seed_;


    //! \brief length of random field
    double periodicity_;

    //! \brief we need Pi
    double pi_;

    //! \brief Do not print to stdout
    bool reduced_output_;

    // deterministic value to use for parameter continuation
    double det_value_paracont_;

    // threshold for fourier series truncation
    int k_trunk_threshold_;

  };  // class RandomFieldFourier
}  // namespace UQ
#endif  /* RANDOMFIELDFOURIER_H_ */
#endif  // HAVE_FFTw
