/*----------------------------------------------------------------------*/
/*! \file
 \brief main file containing routines for calculation of scatra element with chemotactic AND
reactive scalars

\level 2

 *----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_CALC_CHEMO_REAC_HPP
#define BACI_SCATRA_ELE_CALC_CHEMO_REAC_HPP

#include "baci_config.hpp"

#include "baci_mat_scatra_chemotaxis.hpp"
#include "baci_mat_scatra_reaction.hpp"
#include "baci_scatra_ele_calc.hpp"
#include "baci_scatra_ele_calc_advanced_reaction.hpp"
#include "baci_scatra_ele_calc_chemo.hpp"

BACI_NAMESPACE_OPEN


namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
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

BACI_NAMESPACE_CLOSE

#endif  // SCATRA_ELE_CALC_CHEMO_REAC_H
