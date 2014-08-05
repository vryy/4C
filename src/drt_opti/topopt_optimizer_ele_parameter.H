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

  const double Density() {return dens_;}
  const double Viscosity() {return visc_;}

  const double MinPoro() {return min_poro_;}
  const double MaxPoro() {return max_poro_;}
  const double SmearFac() {return smear_fac_;}

  const bool ObjDissipationTerm() {return dissipation_;}
  const bool ObjPresDropTerm() {return pressure_drop_;}
  const double ObjDissipationFac() {return dissipation_fac_;}
  const double ObjPresDropFac() {return pressure_drop_fac_;}

  const bool IsStationary() {return is_stationary_;}
  const INPAR::FLUID::TimeIntegrationScheme TimeIntScheme() {return timealgo_;}
  const bool SUPG() {return supg_;}
  const bool PSPG() {return pspg_;}
  const INPAR::FLUID::TauType TauType() {return whichtau_;}
  const double Dt() {return dt_;}
  const int NumTimesteps() {return max_timesteps_;}
  const double ThetaObj() {return theta_obj_;}
  const double Theta() {return theta_;}

  double VolBd() {return vol_bd_;}

  INPAR::TOPOPT::DensityField DensType() {return dens_type_;}

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
  //! flag to (de-)active supg-stabilization
  bool supg_;
  //! flag to (de-)active pspg-stabilization
  bool pspg_;
  //! Flag to define tau
  INPAR::FLUID::TauType whichtau_;
  //! time-step length
  double dt_;
  /// maximal number of time steps
  int max_timesteps_;
  /// theta in one-step-theta scheme for objective function
  double theta_obj_;
  /// theta in one-step-theta scheme
  double theta_;

  //@}

  double vol_bd_; /// maximal percentage of fluid in background domain

  INPAR::TOPOPT::DensityField dens_type_; /// type of optimization field

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
