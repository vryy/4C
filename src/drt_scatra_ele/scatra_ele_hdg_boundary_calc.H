/*--------------------------------------------------------------------------*/
/*! \file

\brief Routines for ScaTraHDG boundary elements

\level 3

*/
/*--------------------------------------------------------------------------*/


#ifndef SCATRA_ELE_HDG_BOUNDARY_CALC_H
#define SCATRA_ELE_HDG_BOUNDARY_CALC_H

#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


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
    class ScaTraHDGBoundaryImplInterface : DRT::SingletonDestruction
    {
     public:
      //! Empty constructor
      ScaTraHDGBoundaryImplInterface() {}
      //! Empty destructor
      virtual ~ScaTraHDGBoundaryImplInterface() {}

      //! Evaluate a Neumann boundary condition
      /*!
        This class does not provide a definition for this function, it
        must be defined in ScaTraHDGBoundaryImpl.
       */
      virtual int EvaluateNeumann(DRT::ELEMENTS::ScaTraHDGBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseVector& elevec1_epetra) = 0;

      //! Internal implementation class for ScaTraHDGBoundary elements
      static ScaTraHDGBoundaryImplInterface* Impl(const DRT::Element* ele);

    };  // class ScaTraHDGBoundaryImplInterface


    template <DRT::Element::DiscretizationType distype>
    class ScaTraHDGBoundaryImpl : public ScaTraHDGBoundaryImplInterface
    {
     public:
      //! Singleton access method
      static ScaTraHDGBoundaryImpl<distype>* Instance(bool create = true);

      //! called upon destruction
      virtual void Done();

      //! Constructor
      ScaTraHDGBoundaryImpl();

      //! number of element nodes
      static const int bdrynen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions of the ScaTraHDGBoundary element
      static const int bdrynsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! number of space dimensions of the parent element
      static const int nsd_ = bdrynsd_ + 1;

      //! Evaluate a Neumann boundary condition
      virtual int EvaluateNeumann(DRT::ELEMENTS::ScaTraHDGBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseVector& elevec1_epetra);

     private:
      //! node coordinates for boundary element
      LINALG::Matrix<nsd_, bdrynen_> xyze_;
      //! coordinates of current integration point in reference coordinates
      LINALG::Matrix<bdrynsd_, 1> xsi_;
      //! array for shape functions for boundary element
      LINALG::Matrix<bdrynen_, 1> funct_;
      //! array for shape function derivatives for boundary element
      LINALG::Matrix<bdrynsd_, bdrynen_> deriv_;
      //! normal vector pointing out of the domain
      LINALG::Matrix<nsd_, 1> unitnormal_;
      //! velocity vector at integration point
      LINALG::Matrix<nsd_, 1> velint_;
      //! infinitesimal area element drs
      double drs_;
      //! integration factor
      double fac_;

    };  // class ScaTraHDGBoundaryImpl

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* SCATRA_ELE_HDG_BOUNDARY_CALC_H */
