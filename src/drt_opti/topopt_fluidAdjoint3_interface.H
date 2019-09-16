/*---------------------------------------------------------------------*/
/*! \file

\brief Interface for element level of adjoint fluid equations

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef TOPOPT_FLUIDADJOINT3_INTERFACE_H_
#define TOPOPT_FLUIDADJOINT3_INTERFACE_H_


#include "../drt_lib/drt_element.H"


namespace GEO
{
  namespace CUT
  {
    class BoundaryCell;
  }
}  // namespace GEO

namespace DRT
{
  class Discretization;

  namespace UTILS
  {
    class GaussIntegration;
  }

  namespace ELEMENTS
  {
    class Fluid;

    /// Interface base class for FluidAdjoint3Impl
    /*!
      This class exists to provide a common interface for all template
      versions of FluidAdjoint3Impl. The only function
      this class actually defines is Impl, which returns a pointer to
      the appropriate version of FluidAdjoint3Impl.

      \author winklmaier
     */
    class FluidAdjoint3ImplInterface : public DRT::SingletonDestruction
    {
     public:
      /// Empty constructor
      FluidAdjoint3ImplInterface() {}

      /// Evaluate the element
      /*!
        This class does not provide a definition for this function, it
        must be defined in FluidImpl.
       */
      virtual int Evaluate(DRT::ELEMENTS::Fluid* ele, DRT::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
          Epetra_SerialDenseVector& elevec1_epetra) = 0;

      /// Internal implementation class for fluid element
      static FluidAdjoint3ImplInterface* Impl(DRT::Element::DiscretizationType distype);
    };
  }  // namespace ELEMENTS
}  // namespace DRT

#endif /* TOPOPT_FLUIDADJOINT3_INTERFACE_H_ */
