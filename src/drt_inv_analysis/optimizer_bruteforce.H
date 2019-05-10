/*----------------------------------------------------------------------*/
/*!
\brief Brute force sampling of the parameter space

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#ifndef OPTIMIZER_BRUTEFORCE_H_
#define OPTIMIZER_BRUTEFORCE_H_

#include "optimizer_base.H"  // base class

#include "Teuchos_ParameterList.hpp"


namespace INVANA
{
  /*! \class OptimizerBruteForce
   * \brief BruteForce parameter space sampling
   *
   * A nice to have feature to analyse your optimization problem
   * at hand. Of course it is not at all suited to be used in problems
   * with high parametric dimension, e.g. dim>2. At the moment it is
   * restricted to dim=1.
   *
   * Plotting the functional to be minimized at a reduced dimensional
   * problem can help extremely in understanding what is going on.
   *
   * Adjust the hard coded bounds and discretization of the parameter space
   * when using it.
   *
   * \author kehl \date 07/2016
   */
  class OptimizerBruteForce : public OptimizerBase
  {
   public:
    //! constructor
    OptimizerBruteForce(const Teuchos::ParameterList& invp);

    //! destructor
    virtual ~OptimizerBruteForce(){};

    //! brute force sample the parameter space
    virtual void Integrate();

    //! read restart information
    virtual void ReadRestart(int run);

    //! write restart information
    virtual void WriteRestart();

    //! setup algorithm specific stuff
    virtual void Setup();
  };
}  // namespace INVANA



#endif /* OPTIMIZER_BRUTEFORCE_H_ */
