/*!-----------------------------------------------------------------------------
\file   randomvariable.H
\brief Class for generating random variables


 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!---------------------------------------------------------------------------*/
#ifdef HAVE_FFTW
#ifndef RANDOMVARIABLE_H_
#define RANDOMVARIABLE_H_

#include "Teuchos_ParameterList.hpp"
#include <Teuchos_RCP.hpp>

// include random number libraries
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>

// random generator type
typedef boost::mt19937 randnumgen;


namespace UQ
{
  class RandomVariable
  {
   public:
    /*!
    \brief standard constructor

    \param rfp                 (in): parameter list containing the
                                     parameters of pdf etc.

    \param random_variable_id  (in): Id of the random variable in uq problem

    \param seed                (in): some seed for rng

    */
    RandomVariable(const Teuchos::ParameterList& rfp, int random_variable_id, int seed);

    /*!
    \brief destructor

    */
    ~RandomVariable(){};

    /*!
    \brief generates new realization of random variable

    \param seed (in): seed for random number generator

    */
    void CreateNewSample(unsigned int seed);

    /*!
    \brief Evaluate random variable, i.e. return it's value

    \param paracont_parameter (in): parameter continuation parameter [0-1]

    \param writetofile        (in): write some output to file

    \param writetofile        (in): write some output

    */
    double EvalVariable(const double paracont_parameter, const bool writetofile, const bool output);

    /*!
    //! \brief write random field to file

    \param filename        (in): string containing filename including path

    */
    void WriteRandomVariablesToFile(std::string filename, int numrun);

    /*!
    \brief tell me is random variable is gaussian

    */
    bool IsGaussian()
    {
      if (pdf_ == normal)
        return true;
      else
        return false;
    }

    /*!
    \brief do we have a bounded random variable

    If the random variable is bounded EvalVariable() will always
    return values between #lower_bound_ and #upper_bound_

    */
    bool IsBounded() const { return is_bounded_; }

    /*!
    \brief return lower bound of variable

    */
    double LowerBound() const { return lower_bound_; }

    /*!
    \brief return upper bound of variable

    */
    double UpperBound() const { return upper_bound_; }

   private:
    //! \brief possible pdfs
    enum PROBABILITY_DISTRIBUTIONS
    {
      normal,
      lognormal,
      none,
    };

    //! \brief target marginal pdf of random field
    PROBABILITY_DISTRIBUTIONS pdf_;

    //! \brief parameter 1  of probability distribution
    double param_1_;

    //! \brief parameter 2  of probability distribution
    double param_2_;

    //! \brief do we truncate the field above and below a threshold
    bool is_bounded_;

    //! \brief lower cutoff threshold
    double lower_bound_;

    //! \brief upper cutoff threshold
    double upper_bound_;

    //! \brief a random number generator object (boost)
    randnumgen randomnumbergen_;

    //! \brief generates Gaussian random numbers with zero mean and unit standard deviation
    Teuchos::RCP<boost::variate_generator<randnumgen&, boost::normal_distribution<>>> normalgen_;

    //! \brief current value of underlying gaussian random variable
    double value_;

    //! \brief current value of physical parameter
    double physical_value_;

    //! \brief deterministic value to use for parameter continuation
    double det_value_paracont_;

    //! \brief Id of stochastic parameter
    int Id_;

    //! \brief which proc am I on ?
    int myrank_;


  };  // class RandomField
}  // namespace UQ



#endif  /* RANDOMVARIABLE_H_ */
#endif  // FFTW
