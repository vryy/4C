/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_advanced_reaction.H

 \brief main file containing routines for calculation of scatra element with advanced reaction terms


 <pre>
   Maintainer: Moritz Thon
               thon@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_ADVANCED_REACTION_H_
#define SCATRA_ELE_CALC_ADVANCED_REACTION_H_

#include "scatra_ele_calc.H"

namespace DRT
{
namespace ELEMENTS
{

enum reaction_coupling
  {
    simple_multiplicative
  };

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcAdvReac : public virtual ScaTraEleCalc<distype>
{
protected:
    /// (private) protected constructor, since we are a Singleton.
  ScaTraEleCalcAdvReac(const int numdofpernode,const int numscal);

private:

  typedef ScaTraEleCalc<distype> my;

public:
  virtual ~ScaTraEleCalcAdvReac() {}

  /// Singleton access method
  static ScaTraEleCalcAdvReac<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    bool create=true );

  /// called upon destruction
  virtual void Done();

protected:

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
  virtual void MatScaTra(
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

  //! material BioFilm
  virtual void MatBioFilm(
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

  //! material ScatraGrowthScd
  virtual void MatGrowthScd(
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

  virtual void GetRhsInt(
    double&      rhsint,  //!< rhs containing bodyforce at Gauss point
    const double densnp,  //!< density at t_(n+1)
    const int    k        //!< index of current scalar
  );

  //! check if homogenous scatra couplung condition exists and if, read it
  virtual bool IsCoupledAndRead(
      const DRT::Discretization&            ScaTraDiscretization //!<discretisation of the ScaTra field
    );

  virtual double CalcReaCoeff(const int k);

  virtual double CalcReaCoeffFac(
      const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
      const DRT::ELEMENTS::reaction_coupling    couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                                 k                       //!< id of current scalar
  );

  virtual double CalcReaCoeffDerivMatrix(const int k, const int j);

  virtual double CalcReaCoeffDerivFac(
      const std::vector<int>                  stoich,                  //!<stoichometrie of current condition
      const DRT::ELEMENTS::reaction_coupling  couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                               toderive,                //!<concentration to be derived
      const int                               k                       //!< id of current scalar
  );

  virtual double CalcReaBodyForceTerm(const int k);

  virtual double CalcReaBodyForceTermFac(
      const std::vector<int>                      stoich,                 //!<stoichometrie of current condition
      const DRT::ELEMENTS::reaction_coupling      couplingtype            //!<type of coupling the stoichometry coefficients
  );

  virtual double CalcReaBodyForceDerivMatrix(const int k, const int j);

  virtual double CalcReaBodyForceDerivFac(
      const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
      const DRT::ELEMENTS::reaction_coupling    couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                                 toderive                 //!<concentration to be derived
  );

  virtual int Evaluate(DRT::ELEMENTS::Transport*     ele,
                         Teuchos::ParameterList&       params,
                         DRT::Discretization &         discretization,
                         const std::vector<int> &      lm,
                         Epetra_SerialDenseMatrix&     elemat1_epetra,
                         Epetra_SerialDenseMatrix&     elemat2_epetra,
                         Epetra_SerialDenseVector&     elevec1_epetra,
                         Epetra_SerialDenseVector&     elevec2_epetra,
                         Epetra_SerialDenseVector&     elevec3_epetra);

  //! calculation of reactive element matrix
  virtual void CalcMatReact(
    Epetra_SerialDenseMatrix&          emat,        //!< element matrix to be filled
    const int                          k,           //!< index of current scalar
    const double                       timefacfac,  //!< domain-integration factor times time-integration factor
    const double                       timetaufac,  //!< domain-integration factor times time-integration factor times tau
    const double                       taufac,      //!< domain-integration factor times tau
    const double                       densnp,      //!< density at time_(n+1)
    const double                       phinp,       //!< scalar at time_(n+1)
    Teuchos::RCP<ScaTraEleReaManager>  reamanager,  //!< reaction manager
    const LINALG::Matrix<my::nen_,1>&      conv,        //!< convective operator
    const LINALG::Matrix<my::nen_,1>&      sgconv,      //!< subgrid-scale convective operator
    const LINALG::Matrix<my::nen_,1>&      diff         //!< laplace term
  );

  bool iscoupled_;                           /// Flag for homogenous scatra coupling volume condition
  bool isinit_;                              /// Flag if is checked for coupling condition or Material is allready done
  int numcond_;                              /// number of homogenous scatra coupling conditions
  std::vector<std::vector<int> > stoich_;    /// vector containing the stoichometrie
  std::vector<double> reaconst_;             ///constant by which the reaction coefficient is multiplied
  std::vector<DRT::ELEMENTS::reaction_coupling> couplingtype_;///how to calculate the different couplings

}; //end ScaTraEleCalcAdvReac

}

}


#endif /* SCATRA_ELE_CALC_ADVANCED_REACTION_H_ */
