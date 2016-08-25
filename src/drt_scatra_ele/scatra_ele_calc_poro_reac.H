/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_poro_reac.H

 \brief evaluation class containing routines for calculation of scalar transport
        within porous medium including advanced reactions

\level 2

<pre>
  \maintainer Anh-Tu Vuong
              vuong@lnm.mw.tum.de
              http://www.lnm.mw.tum.de
              089 - 289-15264
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
class ScaTraEleCalcPoroReac : virtual public ScaTraEleCalcPoro<distype>, virtual public ScaTraEleCalcAdvReac<distype>
{
protected:
  /// protected constructor, since we are a Singleton.
  ScaTraEleCalcPoroReac(const int numdofpernode,const int numscal,const std::string& disname);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcPoro<distype> poro;
  typedef ScaTraEleCalcAdvReac<distype> advreac;

public:
  virtual ~ScaTraEleCalcPoroReac() {}

  /// Singleton access method
  static ScaTraEleCalcPoroReac<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcPoroReac* delete_me = NULL );

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
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! get the material parameters
  virtual void GetMaterialParams(
    const DRT::Element*   ele,       //!< the element we are dealing with
    std::vector<double>&  densn,     //!< density at t_(n)
    std::vector<double>&  densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>&  densam,    //!< density at t_(n+alpha_M)
    double&               visc,       //!< fluid viscosity
    const int             iquad=-1    //!< id of current gauss point (default = -1)
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

  //! extract element based or nodal values
  //  return extracted values of phinp
  virtual void ExtractElementAndNodeValues(
      DRT::Element*                 ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la
  );

};
}
}


#endif /* SCATRA_ELE_CALC_PORO_H_ */
