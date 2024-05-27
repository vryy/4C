/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation class containing routines for dissolving reactions within
       porous medium for ECM modeling

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_CALC_PORO_REAC_ECM_HPP
#define FOUR_C_SCATRA_ELE_CALC_PORO_REAC_ECM_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc_poro_reac.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  class StructPoroReactionECM;
}


namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype>
    class ScaTraEleCalcPoroReacECM : public ScaTraEleCalcPoroReac<distype>
    {
     private:
      /// private constructor, since we are a Singleton.
      ScaTraEleCalcPoroReacECM(
          const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype> my;
      typedef ScaTraEleCalcPoroReac<distype> pororeac;
      typedef ScaTraEleCalcPoro<distype> poro;
      typedef ScaTraEleCalcAdvReac<distype> advreac;
      using my::nen_;
      using my::nsd_;

     public:
      /// Singleton access method
      static ScaTraEleCalcPoroReacECM<distype>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! get the material parameters
      void get_material_params(const DRT::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,                    //!< density at t_(n)
          std::vector<double>& densnp,                   //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,                   //!< density at t_(n+alpha_M)
          double& visc,                                  //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;

      //! evaluate material
      void materials(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;

      //! compute_struct_chem_potential
      double compute_struct_chem_potential(
          Teuchos::RCP<MAT::StructPoroReactionECM>& structmat, const int gp);
    };

  }  // namespace ELEMENTS
}  // namespace DRT


FOUR_C_NAMESPACE_CLOSE

#endif
