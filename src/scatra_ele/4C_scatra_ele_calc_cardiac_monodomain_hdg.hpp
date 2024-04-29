/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of HDG cardiac monodomain element

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HDG_HPP
#define FOUR_C_SCATRA_ELE_CALC_CARDIAC_MONODOMAIN_HDG_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_polynomial.hpp"
#include "4C_scatra_ele_calc_hdg.hpp"
#include "4C_scatra_ele_hdg.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace ELEMENTS
  {
    /// Scatra HDG element implementation
    template <CORE::FE::CellType distype, int probdim = CORE::FE::dim<distype>>
    class ScaTraEleCalcHDGCardiacMonodomain : public ScaTraEleCalcHDG<distype, probdim>
    {
     protected:
     private:
      /// (private) protected constructor, since we are a Singleton.
      /// this constructor is called from a derived class
      /// -> therefore, it has to be protected instead of private
      ScaTraEleCalcHDGCardiacMonodomain(
          const int numdofpernode, const int numscal, const std::string& disname);

      //    typedef ScaTraEleCalc<distype,probdim> my;
      //    typedef ScaTraEleCalcAniso<distype,probdim> aniso;
      //    typedef ScaTraEleCalcAdvReac<distype,probdim> advreac;

      std::vector<CORE::LINALG::SerialDenseVector> values_mat_gp_all_;
      std::vector<double> gp_mat_alpha_;

     public:
      /// Singleton access method
      static ScaTraEleCalcHDGCardiacMonodomain<distype, probdim>* Instance(const int numdofpernode,
          const int numscal, const std::string& disname, bool create = true);

      //! evaluate the element
      int EvaluateAction(DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, const SCATRA::Action& action,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra)
      {
        return 0;
      };

     protected:
      /*========================================================================*/
      //! @name material and related and related functions
      /*========================================================================*/

      //! evaluate material
      void PrepareMaterials(DRT::Element* ele,  //!< the element we are dealing with
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          Teuchos::RCP<std::vector<CORE::LINALG::SerialDenseMatrix>> difftensor) override;

      //! evaluate material
      virtual void PrepareMaterialsAll(DRT::Element* ele,  //!< the element we are dealing with
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          Teuchos::RCP<std::vector<CORE::LINALG::SerialDenseMatrix>> difftensor);

      //! evaluate material
      virtual void PrepareMaterialsTet(DRT::Element* ele,  //!< the element we are dealing with
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          Teuchos::RCP<std::vector<CORE::LINALG::SerialDenseMatrix>> difftensor);


      //! evaluate material
      void Materials(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          CORE::LINALG::SerialDenseMatrix& difftensor, CORE::LINALG::SerialDenseVector& ivecn,
          CORE::LINALG::SerialDenseVector& ivecnp,
          CORE::LINALG::SerialDenseMatrix& ivecnpderiv) override;

      //! material ScaTra
      void MatMyocard(
          const Teuchos::RCP<const CORE::MAT::Material> material,  //!< pointer to current material
          const int k,                                             //!< id of current scalar
          CORE::LINALG::SerialDenseMatrix& difftensor, CORE::LINALG::SerialDenseVector& ivecn,
          CORE::LINALG::SerialDenseVector& ivecnp, CORE::LINALG::SerialDenseMatrix& ivecnpderiv);

      //! update time dependent material
      void TimeUpdateMaterial(const DRT::Element* ele  //!< the element we are dealing with
          ) override;

      //! get material internal state for output
      void GetMaterialInternalState(const DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization) override;

      //! set material internal state after restart
      void SetMaterialInternalState(const DRT::Element* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization) override;

      //! project material field
      int ProjectMaterialField(const DRT::Element* ele) override;

      //! project material field
      int ProjectMaterialFieldAll(const DRT::Element* ele);

      //! project material field for Tet elements, because quadrature not working for higher order
      //! polynomials with Intrepid
      int ProjectMaterialFieldTet(const DRT::Element* ele);

      /// polynomial space for element interior for various Gauss Points for the evaluation of the
      /// material
      Teuchos::RCP<CORE::FE::PolynomialSpace<probdim>> polySpace_;
    };

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
