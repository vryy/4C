/*--------------------------------------------------------------------------*/
/*! \file

\brief Routines for elemag boundary elements

The routins are copied from the acoustic module and are supposed to be useful in the creation of
boundary elements for the electromagnetic module. The correct implementation is still missing.

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
*/
/*--------------------------------------------------------------------------*/


#ifndef ELEMAG_ELE_BOUNDARY_CALC_H
#define ELEMAG_ELE_BOUNDARY_CALC_H

#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


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
    class ElemagBoundaryImplInterface : DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      ElemagBoundaryImplInterface() {}
      /// Empty destructor
      virtual ~ElemagBoundaryImplInterface() {}

      /// Evaluate a Neumann boundary condition
      /*!
        This class does not provide a definition for this function, it
        must be defined in ElemagBoundaryImpl.
       */
      virtual int EvaluateNeumann(DRT::ELEMENTS::ElemagBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseMatrix* elemat1) = 0;

      /// Evaluate routine for boundary elements inteface
      virtual int Evaluate(DRT::ELEMENTS::ElemagBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra) = 0;

      /// Internal implementation class for ElemagBoundary elements
      static ElemagBoundaryImplInterface* Impl(const DRT::Element* ele);

    };  // class ElemagBoundaryImplInterface


    template <DRT::Element::DiscretizationType distype>
    class ElemagBoundaryImpl : public ElemagBoundaryImplInterface
    {
     public:
      /// Singleton access method
      static ElemagBoundaryImpl<distype>* Instance(bool create = true);

      /// called upon destruction
      virtual void Done();

      /// Constructor
      ElemagBoundaryImpl();

      //! number of element nodes
      static const int bdrynen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      //! number of space dimensions of the ElemagBoundary element
      static const int bdrynsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

      //! number of space dimensions of the parent element
      static const int nsd_ = bdrynsd_ + 1;

      //! Evaluate a Neumann boundary condition
      virtual int EvaluateNeumann(DRT::ELEMENTS::ElemagBoundary* ele,
          Teuchos::ParameterList& params, DRT::Discretization& discretization,
          DRT::Condition& condition, std::vector<int>& lm, Epetra_SerialDenseVector& elevec1_epetra,
          Epetra_SerialDenseMatrix* elemat1);

      /// Evaluate routine for boundary elements
      virtual int Evaluate(DRT::ELEMENTS::ElemagBoundary* ele, Teuchos::ParameterList& params,
          DRT::Discretization& discretization, std::vector<int>& lm,
          Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
          Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
          Epetra_SerialDenseVector& elevec3_epetra);

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

    };  // class ElemagBoundaryImpl

  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* ELEMAG_ELE_BOUNDARY_CALC_H */
