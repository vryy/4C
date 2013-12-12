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

  void SetPlasticParams(Teuchos::ParameterList& params);
  void GetPlasticParams(Teuchos::ParameterList& params);
  bool Converged();
  int NumActivePlasticGP(){return numactive_global_;};
  double DeltaLp_increment_norm() {return lp_increment_norm_global_;};
  double DeltaLp_residual_norm() {return lp_residual_norm_global_;};

private:

  Teuchos::RCP<DRT::Discretization> discret_; // current discretization
  int numactive_global_;// number of active plastic Gauss points on all procs
  bool unconvergedactiveset_; // true, if the active set has not changed in the last step

  double lp_increment_norm_global_; // L2-norm of global vector of plastic flow
  double lp_residual_norm_global_; // L2-norm of global vector of plastic constraint residuals
};// class SsnPlastManager
}//namespace UTILS

#endif // PLASTIC_SSN_MANAGER_H_
