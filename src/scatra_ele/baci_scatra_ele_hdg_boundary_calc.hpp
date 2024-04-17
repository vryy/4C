/*--------------------------------------------------------------------------*/
/*! \file

\brief Routines for ScaTraHDG boundary elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_HDG_BOUNDARY_CALC_HPP
#define FOUR_C_SCATRA_ELE_HDG_BOUNDARY_CALC_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_lib_element.hpp"
#include "baci_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Condition;
  class Discretization;

  namespace ELEMENTS
  {
    class ScaTraHDGBoundary;

    //! Interface base class for ScaTraHDGBoundaryImpl
    /*!
      This class exists to provide a common interface for all template
      versions of ScaTraHDGBoundaryImpl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of ScaTraHDGBoundaryImpl.
     */
    class ScaTraHDGBoundaryImplInterface
    {
     public:
      //! Empty constructor
      ScaTraHDGBoundaryImplInterface() {}
      //! Empty destructor
      virtual ~ScaTraHDGBoundaryImplInterface() = default;
      //! Evaluate a Neumann boundary condition
      /*!
        This class does not provide a definition for this function, it
        must be defined in ScaTraHDGBoundaryImpl.
       */
      virtual int EvaluateNeumann(DRT::ELEMENTS::ScaTraHDGBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra) = 0;

      //! Internal implementation class for ScaTraHDGBoundary elements
      static ScaTraHDGBoundaryImplInterface* Impl(const DRT::Element* ele);

    };  // class ScaTraHDGBoundaryImplInterface


    template <CORE::FE::CellType distype>
    class ScaTraHDGBoundaryImpl : public ScaTraHDGBoundaryImplInterface
    {
     public:
      //! Singleton access method
      static ScaTraHDGBoundaryImpl<distype>* Instance(
          CORE::UTILS::SingletonAction action = CORE::UTILS::SingletonAction::create);

      //! Constructor
      ScaTraHDGBoundaryImpl();

      //! number of element nodes
      static constexpr int bdrynen_ = CORE::FE::num_nodes<distype>;

      //! number of space dimensions of the ScaTraHDGBoundary element
      static constexpr int bdrynsd_ = CORE::FE::dim<distype>;

      //! number of space dimensions of the parent element
      static constexpr int nsd_ = bdrynsd_ + 1;

      //! Evaluate a Neumann boundary condition
      int EvaluateNeumann(DRT::ELEMENTS::ScaTraHDGBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Element::LocationArray& la,
          CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
          CORE::LINALG::SerialDenseVector& elevec1_epetra) override;

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

    };  // class ScaTraHDGBoundaryImpl

  }  // namespace ELEMENTS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
