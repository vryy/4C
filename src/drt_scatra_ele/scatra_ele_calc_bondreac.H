/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_bondreac.H

 \brief main file containing routines for calculation of scatra element with reactive scalars and bond dynamics.

 \level 2

 <pre>
   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289 - 15240
 </pre>
 *----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_BONDREAC_H_
#define SCATRA_ELE_CALC_BONDREAC_H_

#include "scatra_ele_calc.H"
#include "scatra_ele_calc_advanced_reaction.H"

#include "../drt_mat/scatra_reaction_mat.H"
#include "../drt_mat/scatra_bondreac_mat.H"


namespace DRT
{

  // forward declaration
  class ImmersedFieldExchangeManager;

namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype,int probdim=DRT::UTILS::DisTypeToDim<distype>::dim>
class ScaTraEleCalcBondReac : public ScaTraEleCalcAdvReac<distype,probdim>
{
private:
  //! private constructor, since we are a Singleton.
  ScaTraEleCalcBondReac(const int numdofpernode,const int numscal,const std::string& disname);

  typedef ScaTraEleCalc<distype,probdim> my;
  typedef ScaTraEleCalcAdvReac<distype,probdim> advreac;

public:
  virtual ~ScaTraEleCalcBondReac() {}

  //! Singleton access method
  static ScaTraEleCalcBondReac<distype,probdim> * Instance(
      const int numdofpernode,
      const int numscal,
      const std::string& disname,
      const ScaTraEleCalcBondReac* delete_me = NULL );

  //! called upon destruction
  virtual void Done();


protected:

  //! get the material parameters
  virtual void GetMaterialParams(
      const DRT::Element*  ele,       //!< the element we are dealing with
      std::vector<double>& densn,     //!< density at t_(n)
      std::vector<double>& densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
      std::vector<double>& densam,    //!< density at t_(n+alpha_M)
      double&              visc,      //!< fluid viscosity
      const int            iquad      //!< id of current gauss point
  );

  //! set reac. body force, reaction coefficient and derivatives
  virtual void SetBondReactionTerms(
      const int                                 k,           //!< index of current scalar
      const Teuchos::RCP<MAT::MatListBondReacs> matreaclist, //!< index of current scalar
      const double                              violation,   //!< penalty violation at current gp
      const double                              porosity,    //!< porosity of background element
      const double* gpcoord                                  //!< current Gauss-point coordinates
  );


private:

  //! calculate adhesion violation at Gauss point
  double GetAdhesionViolation(
      const DRT::Element*  ele,       //!< the element we are dealing with
      const int            iquad      //!< id of current gauss point
  ) const;

  //! get average distance of integrin to ECM ligand
  double GetPorosityFromBackgroundEle(
      const DRT::Element*  ele,       //!< the element we are dealing with
      const int            iquad      //!< id of current gauss point
  ) const;


}; //end class ScaTraEleCalcBondReac

} //end namespace ELEMENTS

} //end namespace DRT

#endif /* SCATRA_ELE_CALC_BONDREAC_H_ */
