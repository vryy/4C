/*--------------------------------------------------------------------------*/
/*!
\file acou_ele_boundary_calc.H

\brief Routines for acoustic boundary elements

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*--------------------------------------------------------------------------*/


#ifndef ACOU_ELE_BOUNDARY_CALC_H
#define ACOU_ELE_BOUNDARY_CALC_H

#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


namespace DRT
{
class Condition;
class Discretization;

namespace ELEMENTS
{
  class AcouBoundary;

/// Interface base class for AcouBoundaryImpl
/*!
  This class exists to provide a common interface for all template
  versions of AcouBoundaryImpl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of AcouBoundaryImpl.
 */
class AcouBoundaryImplInterface : DRT::SingletonDestruction
{

public:
  /// Empty constructor
  AcouBoundaryImplInterface() {}
  /// Empty destructor
  virtual ~AcouBoundaryImplInterface() {}

  /// Evaluate a Neumann boundary condition
  /*!
    This class does not provide a definition for this function, it
    must be defined in AcouBoundaryImpl.
   */
  virtual int EvaluateNeumann(DRT::ELEMENTS::AcouBoundary*      ele,
      Teuchos::ParameterList&    params,
      DRT::Discretization&       discretization,
      DRT::Condition&            condition,
      std::vector<int>&          lm,
      Epetra_SerialDenseVector&  elevec1_epetra,
      Epetra_SerialDenseMatrix*  elemat1
  ) = 0;

  virtual int Evaluate(
      DRT::ELEMENTS::AcouBoundary*   ele,
      Teuchos::ParameterList&        params,
      DRT::Discretization&           discretization,
      std::vector<int>&              lm,
      Epetra_SerialDenseMatrix&      elemat1_epetra,
      Epetra_SerialDenseMatrix&      elemat2_epetra,
      Epetra_SerialDenseVector&  elevec1_epetra,
      Epetra_SerialDenseVector&  elevec2_epetra,
      Epetra_SerialDenseVector&  elevec3_epetra) = 0;

  /// Internal implementation class for AcouBoundary elements
  static AcouBoundaryImplInterface* Impl(const DRT::Element* ele);

}; // class AcouBoundaryImplInterface


template<DRT::Element::DiscretizationType distype>
class AcouBoundaryImpl: public AcouBoundaryImplInterface
{

public:

  /// Singleton access method
  static AcouBoundaryImpl<distype> * Instance(bool create=true);

  /// called upon destruction
  virtual void Done();

  /// Constructor
  AcouBoundaryImpl();

  //! number of element nodes
  static const int bdrynen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions of the AcouBoundary element
  static const int bdrynsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;

  //! number of space dimensions of the parent element
  static const int nsd_ = bdrynsd_+1;

  //! Evaluate a Neumann boundary condition
  virtual int EvaluateNeumann(DRT::ELEMENTS::AcouBoundary*    ele,
                              Teuchos::ParameterList&    params,
                              DRT::Discretization&       discretization,
                              DRT::Condition&            condition,
                              std::vector<int>&          lm,
                              Epetra_SerialDenseVector&  elevec1_epetra,
                              Epetra_SerialDenseMatrix*  elemat1
                              );

  virtual int Evaluate(
      DRT::ELEMENTS::AcouBoundary*   ele,
      Teuchos::ParameterList&        params,
      DRT::Discretization&           discretization,
      std::vector<int>&              lm,
      Epetra_SerialDenseMatrix&      elemat1_epetra,
      Epetra_SerialDenseMatrix&      elemat2_epetra,
      Epetra_SerialDenseVector&  elevec1_epetra,
      Epetra_SerialDenseVector&  elevec2_epetra,
      Epetra_SerialDenseVector&  elevec3_epetra);

private:
  //! node coordinates for boundary element
  LINALG::Matrix<nsd_,bdrynen_> xyze_;
  //! coordinates of current integration point in reference coordinates
  LINALG::Matrix<bdrynsd_,1> xsi_;
  //! array for shape functions for boundary element
  LINALG::Matrix<bdrynen_,1> funct_;
  //! array for shape function derivatives for boundary element
  LINALG::Matrix<bdrynsd_,bdrynen_> deriv_;
  //! normal vector pointing out of the domain
  LINALG::Matrix<nsd_,1> unitnormal_;
  //! velocity vector at integration point
  LINALG::Matrix<nsd_,1> velint_;
  //! infinitesimal area element drs
  double drs_;
  //! integration factor
  double fac_;

}; // class AcouBoundaryImpl

} // namespace ELEMENTS
} // namespace DRT

#endif /* ACOU_ELE_BOUNDARY_CALC_H */
