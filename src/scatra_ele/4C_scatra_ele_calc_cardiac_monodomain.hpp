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

namespace Discret
{
  namespace ELEMENTS
  {
    template <Core::FE::CellType distype, int probdim>
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
      int evaluate_action(Core::Elements::Element* ele, Teuchos::ParameterList& params,
          Discret::Discretization& discretization, const ScaTra::Action& action,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra,
          Core::LinAlg::SerialDenseVector& elevec2_epetra,
          Core::LinAlg::SerialDenseVector& elevec3_epetra) override;

     protected:
      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! extract element based or nodal values
      void extract_element_and_node_values(Core::Elements::Element* ele,
          Teuchos::ParameterList& params, Discret::Discretization& discretization,
          Core::Elements::Element::LocationArray& la) override;

      //! evaluate material
      void materials(
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
          ) override;

      //! material ScaTra
      void mat_myocard(
          const Teuchos::RCP<const Core::Mat::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          double& densn,                                           //!< density at t_(n)
          double& densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
          double& densam,       //!< density at t_(n+alpha_M)
          double& visc,         //!< fluid viscosity
          const int iquad = -1  //!< id of current gauss point (default = -1)
      );

      //! calculate matrix and rhs for ep
      void sysmat(Core::Elements::Element* ele,       ///< the element whose matrix is calculated
          Core::LinAlg::SerialDenseMatrix& emat,      ///< element matrix to calculate
          Core::LinAlg::SerialDenseVector& erhs,      ///< element rhs to calculate
          Core::LinAlg::SerialDenseVector& subgrdiff  ///< subgrid-diff.-scaling vector
          ) override;
    };

  }  // namespace ELEMENTS

}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
