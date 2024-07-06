/*---------------------------------------------------------------------*/
/*! \file

\brief Methods for randomness

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_UTILS_RANDOM_HPP
#define FOUR_C_UTILS_RANDOM_HPP

#include "4C_config.hpp"

#include <random>

FOUR_C_NAMESPACE_OPEN

namespace Core::UTILS
{
  /*!
  \brief handles random numbers
  */
  class Random
  {
   public:
    /// get a random number
    double uni();

    /// get a vector of random numbers of size count
    void uni(std::vector<double>& randvec, int count);

    /// get a random number
    double normal();

    /// get a vector of random numbers of size count
    void normal(std::vector<double>& randvec, int count);

    /// set the random seed
    void set_rand_seed(const unsigned int seed);

    /// set the range for the uniform rng
    void set_rand_range(const double lower, const double upper);

    /// set the mean and variance for the normal rng
    void set_mean_variance(const double mean, const double var);

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
}  // namespace Core::UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
