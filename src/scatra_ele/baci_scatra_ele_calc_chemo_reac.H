/*----------------------------------------------------------------------*/
/*! \file
 \brief main file containing routines for calculation of scatra element with chemotactic AND
reactive scalars

\level 2

 *----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_CALC_CHEMO_REAC_H
#define BACI_SCATRA_ELE_CALC_CHEMO_REAC_H

#include "baci_mat_scatra_chemotaxis_mat.H"
#include "baci_mat_scatra_reaction_mat.H"
#include "baci_scatra_ele_calc.H"
#include "baci_scatra_ele_calc_advanced_reaction.H"
#include "baci_scatra_ele_calc_chemo.H"


namespace DRT
{
  namespace ELEMENTS
  {
    template <DRT::Element::DiscretizationType distype,
        int probdim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim>
    class ScaTraEleCalcChemoReac : public ScaTraEleCalcChemo<distype, probdim>,
                                   public ScaTraEleCalcAdvReac<distype, probdim>
    {
     private:
      //! private constructor, since we are a Singleton.
      ScaTraEleCalcChemoReac(
          const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype, probdim> my;
      typedef ScaTraEleCalcChemo<distype, probdim> chemo;
      typedef ScaTraEleCalcAdvReac<distype, probdim> advreac;

     public:
      ~ScaTraEleCalcChemoReac() override {}

      //! Singleton access method
      static ScaTraEleCalcChemoReac<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      //! get the material parameters
      void GetMaterialParams(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                  //!< density at t_(n)
          std::vector<double>& densnp,                 //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,                 //!< density at t_(n+alpha_M)
          double& visc,                                //!< fluid viscosity
          const int iquad = -1                         //!< id of current gauss point (default = -1)
          ) override;

    };  // end class ScaTraEleCalcChemoReac

  }  // namespace ELEMENTS

}  // namespace DRT

#endif  // SCATRA_ELE_CALC_CHEMO_REAC_H
