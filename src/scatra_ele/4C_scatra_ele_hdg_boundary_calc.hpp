/*--------------------------------------------------------------------------*/
/*! \file

\brief Routines for ScaTraHDG boundary elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#ifndef FOUR_C_SCATRA_ELE_HDG_BOUNDARY_CALC_HPP
#define FOUR_C_SCATRA_ELE_HDG_BOUNDARY_CALC_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
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
      virtual int evaluate_neumann(Discret::ELEMENTS::ScaTraHDGBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra) = 0;

      //! Internal implementation class for ScaTraHDGBoundary elements
      static ScaTraHDGBoundaryImplInterface* Impl(const Core::Elements::Element* ele);

    };  // class ScaTraHDGBoundaryImplInterface


    template <Core::FE::CellType distype>
    class ScaTraHDGBoundaryImpl : public ScaTraHDGBoundaryImplInterface
    {
     public:
      //! Singleton access method
      static ScaTraHDGBoundaryImpl<distype>* Instance(
          Core::UTILS::SingletonAction action = Core::UTILS::SingletonAction::create);

      //! Constructor
      ScaTraHDGBoundaryImpl();

      //! number of element nodes
      static constexpr int bdrynen_ = Core::FE::num_nodes<distype>;

      //! number of space dimensions of the ScaTraHDGBoundary element
      static constexpr int bdrynsd_ = Core::FE::dim<distype>;

      //! number of space dimensions of the parent element
      static constexpr int nsd_ = bdrynsd_ + 1;

      //! Evaluate a Neumann boundary condition
      int evaluate_neumann(Discret::ELEMENTS::ScaTraHDGBoundary* ele,
          Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
          Core::Elements::Element::LocationArray& la,
          Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
          Core::LinAlg::SerialDenseVector& elevec1_epetra) override;

     private:
      //! node coordinates for boundary element
      Core::LinAlg::Matrix<nsd_, bdrynen_> xyze_;
      //! coordinates of current integration point in reference coordinates
      Core::LinAlg::Matrix<bdrynsd_, 1> xsi_;
      //! array for shape functions for boundary element
      Core::LinAlg::Matrix<bdrynen_, 1> funct_;
      //! array for shape function derivatives for boundary element
      Core::LinAlg::Matrix<bdrynsd_, bdrynen_> deriv_;
      //! normal vector pointing out of the domain
      Core::LinAlg::Matrix<nsd_, 1> unitnormal_;
      //! velocity vector at integration point
      Core::LinAlg::Matrix<nsd_, 1> velint_;
      //! infinitesimal area element drs
      double drs_;
      //! integration factor
      double fac_;

    };  // class ScaTraHDGBoundaryImpl

  }  // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
