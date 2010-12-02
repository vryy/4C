/*----------------------------------------------------------------------*/
/*!
\file fluid3_interface.H

\brief Internal implementation of fluid3 elements

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FLUID3_INTERFACE_H
#define FLUID3_INTERFACE_H

#include "../drt_f3/fluid3.H"

namespace DRT
{
namespace ELEMENTS
{

  /// Interface base class for Fluid3Impl
  /*!
    This class exists to provide a common interface for all template
    versions of Fluid3Impl. The only function
    this class actually defines is Impl, which returns a pointer to
    the appropriate version of Fluid3Impl.
   */
  class Fluid3ImplInterface
  {
  public:

    /// Empty constructor
    Fluid3ImplInterface() {}

    /// Empty destructor
    virtual ~Fluid3ImplInterface() {}

    /// Evaluate the element
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3Impl.
     */
    virtual int Evaluate(DRT::ELEMENTS::Fluid3*    ele,
                         DRT::Discretization & discretization,
                         const std::vector<int> & lm,
                         Teuchos::ParameterList&    params,
                         Teuchos::RCP<MAT::Material> & mat,
                         Epetra_SerialDenseMatrix&  elemat1_epetra,
                         Epetra_SerialDenseMatrix&  elemat2_epetra,
                         Epetra_SerialDenseVector&  elevec1_epetra,
                         Epetra_SerialDenseVector&  elevec2_epetra,
                         Epetra_SerialDenseVector&  elevec3_epetra
      ) = 0;

    /// Integrate shape function
    /*!
      This class does not provide a definition for this function, it
      must be defined in Fluid3Impl.
     */
    virtual int IntegrateShapeFunction(
        DRT::ELEMENTS::Fluid3*    ele,
        DRT::Discretization&      discretization,
        vector<int>&              lm            ,
        Epetra_SerialDenseVector& elevec1       ) = 0;

    /// Internal implementation class for fluid element
    static Fluid3ImplInterface* Impl(DRT::Element::DiscretizationType distype);

  };

}
}

#endif
