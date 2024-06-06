/*---------------------------------------------------------------------*/
/*! \file

\brief Methods for randomness

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_utils_random.hpp"

FOUR_C_NAMESPACE_OPEN

/// get a random number
double Core::UTILS::Random::Uni() { return uni_dist_(rand_engine_); }

/// get a vector of random numbers of size count
void Core::UTILS::Random::Uni(std::vector<double>& randvec, int count)
{
  // resize vector
  randvec.resize(count);

  for (int i = 0; i < count; ++i)
  {
    randvec[i] = uni_dist_(rand_engine_);
  }
}

/// get a random number
double Core::UTILS::Random::Normal() { return norm_dist_(rand_engine_); }

/// get a vector of random numbers of size count
void Core::UTILS::Random::Normal(std::vector<double>& randvec, int count)
{
  // resize vector
  randvec.resize(count);

  for (int i = 0; i < count; ++i)
  {
    randvec[i] = norm_dist_(rand_engine_);
  }
}

/// set the random seed
void Core::UTILS::Random::SetRandSeed(const unsigned int seed) { rand_engine_.seed(seed); }

/// set the range for the uniform rng
void Core::UTILS::Random::SetRandRange(const double lower, const double upper)
{
  std::uniform_real_distribution<double>::param_type parm(lower, upper);
  uni_dist_.param(parm);
}

/// set the mean and variance for the normal rng
void Core::UTILS::Random::SetMeanVariance(const double mean, const double var)
{
  std::normal_distribution<double>::param_type parm(mean, var);
  norm_dist_.param(parm);
}


FOUR_C_NAMESPACE_CLOSE
