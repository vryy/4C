/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.H

\brief evalution of ScaTra elements for ion-transport equation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_ELCH_H
#define SCATRA_ELE_CALC_ELCH_H

#include "scatra_ele_calc.H"

#include "../drt_inpar/inpar_elch.H"

namespace DRT
{
namespace ELEMENTS
{
  // forward declarations
  class ScaTraEleParameterElch;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcElch : public ScaTraEleCalc<distype>
{
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcElch(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;

public:
  virtual ~ScaTraEleCalcElch() {}

  /// Singleton access method
  static ScaTraEleCalcElch<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();

  //! evaluate the element
  int Evaluate(DRT::ELEMENTS::Transport*     ele,
               Teuchos::ParameterList&       params,
               DRT::Discretization &         discretization,
               const std::vector<int> &      lm,
               Epetra_SerialDenseMatrix&     elemat1_epetra,
               Epetra_SerialDenseMatrix&     elemat2_epetra,
               Epetra_SerialDenseVector&     elevec1_epetra,
               Epetra_SerialDenseVector&     elevec2_epetra,
               Epetra_SerialDenseVector&     elevec3_epetra);

  int EvaluateService(DRT::ELEMENTS::Transport*     ele,
                      Teuchos::ParameterList&       params,
                      DRT::Discretization &         discretization,
                      const std::vector<int> &      lm,
                      Epetra_SerialDenseMatrix&     elemat1_epetra,
                      Epetra_SerialDenseMatrix&     elemat2_epetra,
                      Epetra_SerialDenseVector&     elevec1_epetra,
                      Epetra_SerialDenseVector&     elevec2_epetra,
                      Epetra_SerialDenseVector&     elevec3_epetra);

private:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! calculate matrix and rhs. Here the whole thing is hidden.
  void Sysmat(
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    Epetra_SerialDenseVector&             subgrdiff  //!< subgrid-diff.-scaling vector
  );

  /*========================================================================*/
  //! @name material and related and related functions
  /*========================================================================*/

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element* ele,       //!< the element we are dealing with
    double&             densn,     //!< density at t_(n)
    double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    double&             densam,    //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager> diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>  reamanager,   //!< reaction manager
    double&             visc,       //!< fluid viscosity
    const int           iquad=-1    //!< id of current gauss point (default = -1)
  );


  //! evaluate material
  virtual void Materials(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatNewman(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatDiffCond(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatPhase(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               iphase,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  /*========================================================================*/
  //! @name methods for evaluation of individual terms
  /*========================================================================*/


  /*========================================================================*/
  //! @name additional service routines
  /*========================================================================*/

  void CalculateFlux(
  LINALG::Matrix<3,my::nen_>&         flux,
  const DRT::Element*             ele,
  INPAR::ELCH::ElchType           elchtype,
  const INPAR::SCATRA::FluxType   fluxtype,
  const int                       k
  );

  void CalErrorComparedToAnalytSolution(
    const DRT::Element*                   ele,
    const enum INPAR::ELCH::ElchType      elchtype,
    Teuchos::ParameterList&               params,
    Epetra_SerialDenseVector&             errors
  );

  void CalculateConductivity(
    const DRT::Element*               ele,
    const enum INPAR::ELCH::ElchType  elchtype,
    Epetra_SerialDenseVector&         sigma
  );

  void CalculateElectricPotentialField(
    const DRT::Element*               ele,
    const enum INPAR::ELCH::ElchType  elchtype,
    Epetra_SerialDenseMatrix&         emat,
    Epetra_SerialDenseVector&         erhs
  );

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  //! electric potential
  LINALG::Matrix<my::nen_,1> epotnp_;

  //! nodal current density values at t_(n+1) or t_(n+alpha_F)
  LINALG::Matrix<my::nsd_,my::nen_> ecurnp_;

  DRT::ELEMENTS::ScaTraEleParameterElch* elchpara_;

  INPAR::ELCH::DiffCondMat diffcondmat_;

  bool cursolvar_;

  enum INPAR::ELCH::EquPot equpot_;

  //! flag for migration velocity in the stabilization operator (ELCH)
  bool migrationstab_;
  //! flag for migration velocity in formula for stabilization parameter (ELCH)
  bool migrationintau_;
  //! flag for migration velocity in residual of strong form (ELCH)
  bool migrationinresidual_;

  //TODO: BACI_ELE_CLEANING
  //! charge numbers of all ion species (ELCH specific)
  std::vector<double> valence_;

  //! transference numbers
  std::vector<double> trans_;
  double transelim_;
  std::vector<std::vector<double> > transderiv_;

  //! conductivity
  std::vector<double> cond_;
  std::vector<double> condderiv_;

  //! thermodynamic factor
  std::vector<double> therm_;
  std::vector<double> thermderiv_;

  std::vector<double> diffusderiv_;
  std::vector<double> diffuselimderiv_;
  double diffuselim_;
  std::vector<double> eps_;
  std::vector<double> tort_;
  std::vector<double> epstort_;

  /// dilute solution theory (diffusion potential in current equation):
  ///    a          b
  ///   |--|  |----------|
  ///   z_1 + (z_2 - z_1) t_1
  /// ------------------------ (RT/F kappa 1/c_k grad c_k)
  ///      z_1 z_2
  ///     |________|
  ///         c

  double a_;
  double b_;
  double c_;

};
}
}

#endif

