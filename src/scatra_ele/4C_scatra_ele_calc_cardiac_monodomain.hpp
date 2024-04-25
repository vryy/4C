/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_cardiac_monodomain.H

\level 2

 *----------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HPP
#define FOUR_C_SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HPP

#include "4C_config.hpp"

#include "4C_scatra_ele_calc.hpp"
#include "4C_scatra_ele_calc_advanced_reaction.hpp"
#include "4C_scatra_ele_calc_aniso.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    template <CORE::FE::CellType distype, int probdim>
    class ScaTraEleCalcCardiacMonodomain : public ScaTraEleCalcAniso<distype, probdim>,
                                           public ScaTraEleCalcAdvReac<distype, probdim>
    {
     protected:
     private:
      /// (private) protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcCardiacMonodomain(
          const int numdofpernode, const int numscal, const std::string& disname);

      typedef ScaTraEleCalc<distype, probdim> my;
      typedef ScaTraEleCalcAniso<distype, probdim> aniso;
      typedef ScaTraEleCalcAdvReac<distype, probdim> advreac;
      using my::nen_;
      using my::nsd_;
      using my::nsd_ele_;

     public:
      /// Singleton access method
      static ScaTraEleCalcCardiacMonodomain<distype, probdim>* Instance(
          const int numdofpernode, const int numscal, const std::string& disname);

      //! evaluate the element
      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

     protected:
      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! extract element based or nodal values
      void ExtractElementAndNodeValues(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la) override;

      //! evaluate material
      void Materials(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;

      //! material ScaTra
      void MatMyocard(
          const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
          const int k,                                       //!< id of current scalar
          double& densn,                                     //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! calculate matrix and rhs for ep
      void Sysmat(DRT::Element* ele,                  ///< the element whose matrix is calculated
          CORE::LINALG::SerialDenseMatrix& emat,      ///< element matrix to calculate
          CORE::LINALG::SerialDenseVector& erhs,      ///< element rhs to calculate
          CORE::LINALG::SerialDenseVector& subgrdiff  ///< subgrid-diff.-scaling vector
          ) override;
    };

  }  // namespace ELEMENTS

}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
