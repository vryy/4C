/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer_ele_parameter.H

\brief Evaluation of element parameter

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
#include "../drt_inpar/inpar_topopt.H"


namespace DRT
{
namespace ELEMENTS
{

/*!
\brief Evaluation of element parameter

This class just sets most of the element parameter of the topopt-element in a
pre-evaluation step due to performance reasons

\author Martin Winklmaier (winklmaier@lnm.mw.tum.de)
*/class TopOptParam
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

  double Density() {return dens_;}
  double Viscosity() {return visc_;}

  double MinPoro() {return min_poro_;}
  double MaxPoro() {return max_poro_;}
  double SmearFac() {return smear_fac_;}

  bool ObjDissipationTerm() {return dissipation_;}
  bool ObjPresDropTerm() {return pressure_drop_;}
  double ObjDissipationFac() {return dissipation_fac_;}
  double ObjPresDropFac() {return pressure_drop_fac_;}
//
  bool IsStationary() {return is_stationary_;}
  INPAR::FLUID::TimeIntegrationScheme TimeIntScheme() {return timealgo_;}
  INPAR::FLUID::TauType TauType() {return whichtau_;}
  double Dt() {return dt_;}
  int NumTimesteps() {return max_timesteps_;}
  double Theta() {return theta_;}
  double ThetaPre() {return theta_pre_;}
  double ThetaDiv() {return theta_div_;}

  double VolBd() {return vol_bd_;}

  INPAR::TOPOPT::OptiCase OptiCase() {return opti_case_;}

private:

  /// @name flow material parameters

  //! fluid density
  double dens_;
  /// fluid viscosity
  double visc_;

  //@}


  /// @name optimization material parameters

  double min_poro_;
  double max_poro_;
  double smear_fac_;

  //@}


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
  //! Flag to define tau
  INPAR::FLUID::TauType whichtau_;
  //! time-step length
  double dt_;
  /// maximal number of time steps
  int max_timesteps_;
  /// theta in one-step-theta scheme
  double theta_;
  /// theta for pressure terms in one-step-theta scheme
  double theta_pre_;
  /// theta for divergence terms in one-step-theta scheme
  double theta_div_;

  //@}

  double vol_bd_; /// maximal percentage of fluid in background domain

  INPAR::TOPOPT::OptiCase opti_case_; /// type of optimization (e.g. special test cases)

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
