/*!-------------------------------------------------------------------------*\
 * \file meshfree_scatra_impl.H
 *
 * \brief Interface of meshfree scalar transport problems
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#ifndef MESHFREE_SCATRA_CELL_INTERFACE_H
#define MESHFREE_SCATRA_CELL_INTERFACE_H

#include <vector>
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_singletondestruction.H"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

namespace DRT
{
class Discretization;

namespace ELEMENTS
{
class MeshfreeTransport;

/*--------------------------------------------------------------------------*/
/*!
 * \brief Interface base class for ScaTraImpl
 *
 * This class exists to provide a common interface for all template
 * versions of ScaTraImpl. The only function
 * this class actually defines is Impl, which returns a pointer to
 * the appropriate version of ScaTraImpl.
 * Deriving from DRT::SingletonDestruction provides proper registration
 * of singleton at DRT::Problem instance and - as a result -
 * correct memory handling at the end of a BACI run.
 * (correctly delete allocated ScaTraImplInterface instances)
 *
 * \author Keijo Nissen (nissen@lnm.mw.tum.de)
 *
 * \date March, 2012
 */
/*--------------------------------------------------------------------------*/
class MeshfreeScaTraCellInterface : public DRT::SingletonDestruction
{

public:
  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Empty constructor
   *///                                                   (public) nis Jan12
  /*------------------------------------------------------------------------*/
  MeshfreeScaTraCellInterface() {}

  /*========================================================================*/
  //! @name Evaluate methods
  /*========================================================================*/

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Evaluate the cell - virtual function
   *
   * This class does not provide a definition for this function, it must be
   * defined in ScaTraImpl.
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  virtual int Evaluate(
    DRT::ELEMENTS::MeshfreeTransport* ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization &             discretization,
    const std::vector<int> &          lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseMatrix&         elemat2_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    Epetra_SerialDenseVector&         elevec2_epetra,
    Epetra_SerialDenseVector&         elevec3_epetra) = 0;

  virtual int EvaluateService(
    DRT::ELEMENTS::MeshfreeTransport* ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization &             discretization,
    const std::vector<int> &          lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseMatrix&         elemat2_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    Epetra_SerialDenseVector&         elevec2_epetra,
    Epetra_SerialDenseVector&         elevec3_epetra) = 0;

};

}
}
#endif // MESHFREE_SCATRA_CELL_INTERFACE_H
