/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer_ele_parameter.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#ifndef TOPOPT_OPTIMIZER_ELE_PARAMETER_H_
#define TOPOPT_OPTIMIZER_ELE_PARAMETER_H_


#include "../drt_inpar/inpar_fluid.H"


namespace DRT
{
namespace ELEMENTS
{

/// Evaluation of general parameters
class TopOptParam
{
public:

  //! Singleton access method
  static Teuchos::RCP<TopOptParam> Instance();

  //! general fluid parameter are set
  void SetGeneralOptimizationParameter(
      Teuchos::ParameterList& params    //> parameter list
  );

  //! print parameter to screen
  void PrintAdjointParameter() const;

  bool ObjDissipationTerm() {return dissipation_;}
  bool ObjPresDropTerm() {return pressure_drop_;}
  double ObjDissipationFac() {return dissipation_fac_;}
  double ObjPresDropFac() {return pressure_drop_fac_;}
//
  bool IsStationary() {return is_stationary_;}
  INPAR::FLUID::TimeIntegrationScheme TimeIntScheme() {return timealgo_;}
  double Dt() {return dt_;}

private:

  /// @name objective parameters
  // set if objective contains dissipation
  bool dissipation_;
  // set if objective contains pressure drop
  bool pressure_drop_;
  // set objective's dissipation factor
  double dissipation_fac_;
  // set objective's pressure drop factor
  double pressure_drop_fac_;

  //@}

  /// @name parameters for the time integration
  //! flag to (de)activate stationary formulation
  bool is_stationary_;
  //! time algorithm
  INPAR::FLUID::TimeIntegrationScheme timealgo_;
  //! time-step length
  double dt_;
  /// maximal number of time steps
  int max_timesteps_;

  //@}

  // private constructor
  TopOptParam();

  /// no copying
  TopOptParam(const TopOptParam&);

  /// no copying
  TopOptParam& operator=(const TopOptParam&);

  /// the single instance
  static Teuchos::RCP<TopOptParam> instance_;
};

}
}

#endif /* TOPOPT_OPTIMIZER_ELE_PARAMETER_H_ */
