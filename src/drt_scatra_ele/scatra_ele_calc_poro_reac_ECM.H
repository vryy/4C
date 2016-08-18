/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_poro_reac_ECM.H

\brief scatra_ele_calc_poro_reac_ECM.H

\level 2

<pre>
  \maintainer Anh-Tu Vuong
              vuong@lnm.mw.tum.de
              http://www.lnm.mw.tum.de
              089 - 289-15251
</pre>
 *----------------------------------------------------------------------*/


#ifndef SCATRA_ELE_CALC_PORO_REAC_ECM_H_
#define SCATRA_ELE_CALC_PORO_REAC_ECM_H_

#include "scatra_ele_calc_poro_reac.H"

namespace MAT
{
class StructPoroReactionECM;
}


namespace DRT
{
namespace ELEMENTS
{

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcPoroReacECM : public ScaTraEleCalcPoroReac<distype>
{
private:
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcPoroReacECM(const int numdofpernode,const int numscal,const std::string& disname);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcPoroReac<distype> pororeac;
  typedef ScaTraEleCalcPoro<distype> poro;
  typedef ScaTraEleCalcAdvReac<distype> advreac;

public:
  virtual ~ScaTraEleCalcPoroReacECM() {}

  /// Singleton access method
  static ScaTraEleCalcPoroReacECM<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcPoroReacECM* delete_me = NULL );

  /// called upon destruction
  virtual void Done();

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element* ele,       //!< the element we are dealing with
    double&             densn,     //!< density at t_(n)
    double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    double&             densam,    //!< density at t_(n+alpha_M)
    double&             visc,       //!< fluid viscosity
    const int           iquad=-1    //!< id of current gauss point (default = -1)
  );

  //! Get numcond, stoich list, reaction coefficient, couplingtpye from material
  virtual void GetAdvancedReactionCoefficients(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      const int           iquad
  );

  //! evaluate material
  virtual void Materials(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! ComputeStructChemPotential
  double ComputeStructChemPotential(
      Teuchos::RCP<MAT::StructPoroReactionECM>& structmat,
      const int gp);
};

/// Scatra reaction manager for Advanced_Reaction
/*!
  This class keeps all Advanced_Reaction-specific reaction parameter needed for the evaluation of an element. The
  ScaTraEleReaManagerAdvReac is derived from the standard ScaTraEleReaManager.
*/
class ScaTraEleReaManagerPoroReacECM : public ScaTraEleReaManagerAdvReac
{
public:

  ScaTraEleReaManagerPoroReacECM(int numscal): ScaTraEleReaManagerAdvReac(numscal)
    {return;}

  virtual ~ScaTraEleReaManagerPoroReacECM(){};
};

}
}


#endif /* SCATRA_ELE_CALC_PORO_REAC_ECM_H_ */
