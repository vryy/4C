/*---------------------------------------------------------------------*/
/*! \file

\brief Methods for randomness

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_UTILS_RANDOM_HPP
#define FOUR_C_UTILS_RANDOM_HPP

#include "baci_config.hpp"

#include <random>

BACI_NAMESPACE_OPEN

namespace CORE::UTILS
{
  /*!
  \brief handles random numbers
  */
  class Random
  {
   public:
    /// get a random number
    double Uni();

    /// get a vector of random numbers of size count
    void Uni(std::vector<double>& randvec, int count);

    /// get a random number
    double Normal();

    /// get a vector of random numbers of size count
    void Normal(std::vector<double>& randvec, int count);

    /// set the random seed
    void SetRandSeed(const unsigned int seed);

    /// set the range for the uniform rng
    void SetRandRange(const double lower, const double upper);

    /// set the mean and variance for the normal rng
    void SetMeanVariance(const double mean, const double var);

   private:
    /// @name Random number generation
    /// @{
    /// random number engine
    std::mt19937 rand_engine_{0};

    /// uniform random number distribution between -1.0 and 1.0
    std::uniform_real_distribution<double> uni_dist_{-1.0, 1.0};

    /// unit normal random number distribution
    std::normal_distribution<double> norm_dist_{};
    //@}
  };
}  // namespace CORE::UTILS

BACI_NAMESPACE_CLOSE

#endif
