/*!----------------------------------------------------------------------
\file plastic_ssn_manager.H

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#ifndef PLASTIC_SSN_MANAGER_H_
#define PLASTIC_SSN_MANAGER_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace UTILS
{
/*!
\brief Class controlling semi-smooth Newton strategy for plasticity
*/
class PlastSsnManager
{
public:

  /*!
    \brief Standard Constructor
  */
  PlastSsnManager(Teuchos::RCP<DRT::Discretization> discret);

  /*!
    \brief Destructor
  */
  ~PlastSsnManager() {};

  //@}

  //! @name Evaluation methods

  /*!
  \brief Read and check plasticity input parameters

  All specified plasticity-related input parameters are read from the
  DRT::Problem::Instance() and stored into a local variable of
  type Teuchos::ParameterList. Invalid parameter combinations are
  sorted out and throw a dserror.

  */
  void ReadAndCheckInput();
  
  /*
  \brief set plastic parameters before element evaluation
  */
  void SetPlasticParams(Teuchos::ParameterList& params);

  /*
   \brief get plastic parameters from elements after evaluation
   */
  void GetPlasticParams(Teuchos::ParameterList& params);

  /*
  \brief check active set for convergence
  */
  bool ActiveSetConverged();

  /*
  \brief check constraint residual for convergence
  */
  bool ConstraintConverged();

  /*
  \brief check plastic flow increment for convergence
  */
  bool IncrementConverged();

  /*
  \brief return total number of active Gauss points
  */
  int NumActivePlasticGP(){return numactive_global_;};

  /*
  \brief return L2 norm of plastic flow increment
  */
  double DeltaLp_increment_norm() {return lp_increment_norm_global_;};

  /*
  \brief return L2 norm of constraint equations
  */
  double DeltaLp_residual_norm() {return lp_residual_norm_global_;};

  /*
  \brief return parameter list of plastic algorithmic parameters
  */
  Teuchos::RCP<Teuchos::ParameterList> Params() {return plparams_; };

private:

  Teuchos::RCP<DRT::Discretization> discret_; // current discretization
  Teuchos::RCP<Teuchos::ParameterList> plparams_; // plasticity algorithm input parameters

  int numactive_global_;// number of active plastic Gauss points on all procs
  bool unconvergedactiveset_; // true, if the active set has not changed in the last step

  double lp_increment_norm_global_; // L2-norm of global vector of plastic flow
  double lp_residual_norm_global_; // L2-norm of global vector of plastic constraint residuals
};// class SsnPlastManager
}//namespace UTILS

#endif // PLASTIC_SSN_MANAGER_H_
