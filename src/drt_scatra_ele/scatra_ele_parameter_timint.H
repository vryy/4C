/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_parameter_timint.H

\brief Setting of time-integration parameters in scatra element evaluation

This class provides the scatra element parameter for the time integration,
which are unique and equal for every scatra in the problem. Time integration
with different parameters in more than one scatra field is not yet supported.

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_PARAMETER_TIMINT_H
#define SCATRA_ELE_PARAMETER_TIMINT_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_singletondestruction.H"

namespace DRT
{
namespace ELEMENTS
{

/// Evaluation of general parameters (constant over time)
class ScaTraEleParameterTimInt : public DRT::SingletonDestruction
{
private:

  // private constructor
  ScaTraEleParameterTimInt();

public:

  //! Singleton access method
  static ScaTraEleParameterTimInt* Instance(const std::string& disname, bool create=true);

  //! Singleton destruction method called from global problem
  void Done();

  //! dtor
  ~ScaTraEleParameterTimInt() {}

  //! time parameter are set
  void SetElementTimeParameter(Teuchos::ParameterList& params    //> parameter list
      );

  //! print parameter to screen
  void PrintScaTraTimeParameter();

  bool IsGenAlpha(){return is_genalpha_;};
  bool IsStationary(){return is_stationary_;};
  bool IsIncremental(){return is_incremental_;};
  double Time(){return time_;};
  double Dt(){return dt_;};
  double TimeFac(){return timefac_;};
  double TimeFacRhs(){return timefacrhs_;};
  double TimeFacRhsTau(){return timefacrhstau_;};
  double AlphaF(){return alphaF_;};


private:

  bool is_genalpha_;
  bool is_stationary_;
  bool is_incremental_;
  double time_;
  double dt_;
  double timefac_;
  double timefacrhs_;
  double timefacrhstau_;
//  theta_(0.0),
//  omtheta_(0.0),
//  gamma_(0.0),
  double alphaF_;
//  alphaM_(0.0),
//  afgdt_(1.0),
//  timefacrhs_(1.0),
//  timefacpre_(1.0)


}; // class ScaTraEleParameterTimInt

} // namespace ELEMENTS
} // namespace DRT

#endif
