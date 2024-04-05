/*--------------------------------------------------------------------------*/
/*! \file

\brief Routines for elemag boundary elements

The routines are used in the creation of boundary elements for the electromagnetic module. The
correct implementation is still missing.

\level 2

*/
/*--------------------------------------------------------------------------*/


#ifndef FOUR_C_ELEMAG_ELE_BOUNDARY_CALC_HPP
#define FOUR_C_ELEMAG_ELE_BOUNDARY_CALC_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_lib_element.hpp"
#include "baci_utils_singleton_owner.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Condition;
  class Discretization;

  namespace ELEMENTS
  {
    class ElemagBoundary;

    /// Interface base class for ElemagBoundaryImpl
    /*!
      This class exists to provide a common interface for all template
      versions of ElemagBoundaryImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of ElemagBoundaryImpl.
     */
    class ElemagBoundaryImplInterface
    {
     public:
      /// Empty constructor
      ElemagBoundaryImplInterface() {}
      /// Empty destructor
      virtual ~ElemagBoundaryImplInterface() = default;
      /// Evaluate a Neumann boundary condition
      /*!
        This class does not provide a definition for this function, it
        must be defined in ElemagBoundaryImpl.
       */
      virtual int EvaluateNeumann(DRT::ELEMENTS::ElemagBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseMatrix* elemat1) = 0;

      /// Evaluate routine for boundary elements inteface
      virtual int Evaluate(DRT::ELEMENTS::ElemagBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) = 0;

      /// Internal implementation class for ElemagBoundary elements
      static ElemagBoundaryImplInterface* Impl(const DRT::Element* ele);

    };  // class ElemagBoundaryImplInterface


    template <CORE::FE::CellType distype>
    class ElemagBoundaryImpl : public ElemagBoundaryImplInterface
    {
     public:
      /// Singleton access method
      static ElemagBoundaryImpl<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      /// Constructor
      ElemagBoundaryImpl();

      //! number of element nodes
      static constexpr int bdrynen_ = CORE::FE::num_nodes<distype>;

      //! number of space dimensions of the ElemagBoundary element
      static constexpr int bdrynsd_ = CORE::FE::dim<distype>;

      //! number of space dimensions of the parent element
      static constexpr int nsd_ = bdrynsd_ + 1;

      //! Evaluate a Neumann boundary condition
      int EvaluateNeumann(DRT::ELEMENTS::ElemagBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseMatrix* elemat1) override;

      /// Evaluate routine for boundary elements
      int Evaluate(DRT::ELEMENTS::ElemagBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra,
          CORE::LINALG::SerialDenseVector& elevec2_epetra,
          CORE::LINALG::SerialDenseVector& elevec3_epetra) override;

     private:
      //! node coordinates for boundary element
      CORE::LINALG::Matrix<nsd_, bdrynen_> xyze_;
      //! coordinates of current integration point in reference coordinates
      CORE::LINALG::Matrix<bdrynsd_, 1> xsi_;
      //! array for shape functions for boundary element
      CORE::LINALG::Matrix<bdrynen_, 1> funct_;
      //! array for shape function derivatives for boundary element
      CORE::LINALG::Matrix<bdrynsd_, bdrynen_> deriv_;
      //! normal vector pointing out of the domain
      CORE::LINALG::Matrix<nsd_, 1> unitnormal_;
      //! velocity vector at integration point
      CORE::LINALG::Matrix<nsd_, 1> velint_;
      //! infinitesimal area element drs
      double drs_;
      //! integration factor
      double fac_;

    };  // class ElemagBoundaryImpl

  }  // namespace ELEMENTS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // ELEMAG_ELE_BOUNDARY_CALC_H
