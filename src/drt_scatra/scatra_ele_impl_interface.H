/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_impl_interface.H

\brief General layout for internal implementation classes for
       scalar transport element routines

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_IMPL_INTERFACE_H
#define SCATRA_ELE_IMPL_INTERFACE_H

#include "../drt_lib/drt_singletondestruction.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_inpar/inpar_scatra.H"

namespace DRT
{
namespace ELEMENTS
{

//! Interface base class for ScaTraImpl
/*!
  This class exists to provide a common interface for all template
  versions of ScaTraImpl. The only function
  this class actually defines is Impl, which returns a pointer to
  the appropriate version of ScaTraImpl.
  Deriving from DRT::SingletonDestruction provides proper registration
  of singleton at DRT::Problem instance and - as a result -
  correct memory handling at the end of a BACI run.
  (correctly delete allocated ScaTraImplInterface instances)
 */
class ScaTraImplInterface : public DRT::SingletonDestruction
{
public:
  //! Empty constructor
  ScaTraImplInterface() {}
  //! Empty destructor
  virtual ~ScaTraImplInterface() {}
  //! Evaluate the element
  /*!
    This class does not provide a definition for this function, it
    must be defined in ScaTraImpl.
   */
  virtual int Evaluate(DRT::Element*              ele,
                       Teuchos::ParameterList&    params,
                       DRT::Discretization&       discretization,
                       std::vector<int>&          lm,
                       Epetra_SerialDenseMatrix&  elemat1_epetra,
                       Epetra_SerialDenseMatrix&  elemat2_epetra,
                       Epetra_SerialDenseVector&  elevec1_epetra,
                       Epetra_SerialDenseVector&  elevec2_epetra,
                       Epetra_SerialDenseVector&  elevec3_epetra
                       ) = 0;

  //! return appropriate internal implementation class (acts as a simple factory)
  static ScaTraImplInterface* Impl(
      const DRT::Element* ele,
      const enum INPAR::SCATRA::ScaTraType scatratype,
      const bool reinit
      );

}; // class ScaTraImplInterface

} // namespace ELEMENTS
} // namespace DRT


#endif // SCATRA_ELE_IMPL_INTERFACE_H
