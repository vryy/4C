/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro_reac.H

 \brief

 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_PORO_REAC_H_
#define SCATRA_ELE_CALC_PORO_REAC_H_

#include "scatra_ele_calc.H"
#include "scatra_ele_calc_poro.H"
#include "scatra_ele_calc_advanced_reaction.H"

namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcPoroReac : public ScaTraEleCalcPoro<distype>, public ScaTraEleCalcAdvReac<distype>
{
private:
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcPoroReac(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcPoro<distype> poro;
  typedef ScaTraEleCalcAdvReac<distype> advreac;

public:
  virtual ~ScaTraEleCalcPoroReac() {}

  /// Singleton access method
  static ScaTraEleCalcPoroReac<distype> * Instance(
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

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element* ele,       //!< the element we are dealing with
    double&             densn,     //!< density at t_(n)
    double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    double&             densam,    //!< density at t_(n+alpha_M)
    double&             visc,       //!< fluid viscosity
    const int           iquad=-1    //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatScaTra(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material ScaTra
  virtual void MatPoroECM(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

private:

  void SetReactionTermsMatScatra(const int                               k,
                                const double                            porosity);

  //!Calculate K(c)
  double CalcReaCoeff(
      const int                        k,                      //!< id of current scalar
      const double                     porosity                //!< current porosity
      );


  double CalcReaCoeffFac(
      const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
      const MAT::PAR::reaction_coupling         couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                                 k,                       //!< id of current scalar
      const double                              porosity                 //!< current porosity
  );


  //! calculate \frac{partial}{\partial c} K(c)
  double CalcReaCoeffDerivMatrix(
      const int                 k,                  //!< id of current scalar
      const int                 j,                  //!< concentration to be derived to
      const double              porosity            //!< current porosity
      );


  double CalcReaCoeffDerivFac(
      const std::vector<int>                  stoich,                  //!<stoichometrie of current condition
      const MAT::PAR::reaction_coupling       couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                               toderive,                //!<concentration to be derived to
      const int                               k,                       //!< id of current scalar
      const double                            porosity                //!< current porosity
  );


  //! calculate f(c)
  double CalcReaBodyForceTerm(
      const int                                 k,                      //!< id of current scalar
      const double                              porosity                //!< current porosity
);


  double CalcReaBodyForceTermFac(
      const std::vector<int>                      stoich,                 //!<stoichometrie of current condition
      const MAT::PAR::reaction_coupling           couplingtype,           //!<type of coupling the stoichometry coefficients
      const double                                porosity                //!< current porosity
  );


  //! calculate \frac{partial}{\partial c} f(c)
  double CalcReaBodyForceDerivMatrix(
      const int                 k,                  //!< id of current scalar
      const int                 j,                  //!< concentration to be derived to
      const double              porosity                //!< current porosity
);


  virtual double CalcReaBodyForceDerivFac(
      const std::vector<int>                    stoich,                  //!<stoichometrie of current condition
      const MAT::PAR::reaction_coupling         couplingtype,            //!<type of coupling the stoichometry coefficients
      const int                                 toderive,                //!<concentration to be derived to
      const double                              porosity                //!< current porosity
  );

};
}
}


#endif /* SCATRA_ELE_CALC_PORO_H_ */
