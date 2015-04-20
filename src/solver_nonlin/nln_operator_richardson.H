/*----------------------------------------------------------------------------*/
/*!
 \file nln_operator_richardson.H

 \brief Richardson iterations with relaxation and preconditioning

 <pre>
 Maintainer: Matthias Mayr
             mayr@mhpc.mw.tum.de
             089 - 289-10362
 </pre>
 */

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_RICHARDSON_H
#define NLN_OPERATOR_RICHARDSON_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "nln_operator_base.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */

class Epetra_Comm;
class Epetra_MultiVector;

namespace LINALG
{
namespace SOLVER
{
class PreconditionerType;
} // namespace SOLVER
} // namespace LINALG

namespace NLNSOL
{
class LineSearchBase;
} // namespace NLNSOL

namespace Teuchos
{
class ParameterList;
} // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
/*! \brief Richardson iteration with relaxation and preconditioning
 *
 *  This algorithm implements a Richardson iteration with IFPACK
 *  preconditioning. We allow for a user defined relaxation parameter.
 *
 *  As IFPACK preconditioners, we allow for
 *  <ul>
 *  <li> Jacobi </li>
 *  <li> Symmetric Gauss-Seidel </li>
 *  <li> Chebyshev </li>
 *  </ul>
 *  In case of relaxation-type preconditioning, set the number of sweeps to 1
 *  for the IFPACK configuration. Sweeps/iterations as well as the damping are
 *  performed by the Richardson algorithm.
 *
 *  \sa NlnOperatorBase
 *
 *  \author mayr.mt \date 04/2015
 */
class NlnOperatorRichardson : public NlnOperatorBase
{
public:

  //! @name Construction
  //@{

  //! Constructor (empty)
  NlnOperatorRichardson();

  //! Destructor
  virtual ~NlnOperatorRichardson()
  {
  }

  //@}

  //! @name Setup
  //@{

  //! Setup of the preconditioner object
  virtual void Setup();

  /*! \brief Rebuild preconditioner
   *
   *  No action required.
   *
   *  \author mayr.mt \date 04/2015
   */
  virtual void RebuildPrec(){ /* no action required */ };

  //@}

  //! @name Mathematical functions
  //@{

  /*! \brief Apply the preconditioner
   *
   *
   *
   *  \return Integer error code, set to 0 if successful.
   *
   *  \author mayr.mt \date 04/2015
   */
  virtual int ApplyInverse(const Epetra_MultiVector& f, ///< residual in mechanical sense
      Epetra_MultiVector& x ///< solution
      ) const;

  //! @name Attribute access functions
  //@{

  //! Returns a character string describing the operator
  virtual const char* Label() const { return "NLNSOL::NlnOperatorRichardson"; }

  //! Returns a character string describing the operator briefly
  virtual const char* LabelShort() const { return "Richardson"; }

  //@}

protected:

private:

  //! @name Search Direction
  //@{

  //! Compute line search parameter (=step length)
  virtual void ComputeStepLength(
      const Epetra_MultiVector& x, ///< current solution
      const Epetra_MultiVector& f, ///< residual at current solution
      const Epetra_MultiVector& inc, ///< search direction
      const double fnorm2, ///< current L2-norm of residual
      double& lsparam, ///< line search parameter (to be determined)
      bool& suffdecr ///< achieved sufficient decrease? (to be determined)
      ) const;

  //@}

  /*! \brief Create the IFPACK parameter list
   *
   *  The IFPACK parameter list is read from a SOLVER block from the Baci input
   *  file.
   *
   *  Some additional checks are performed. Especially, we have to ensure that
   *  the IFPACK preconditioner is applied only once when calling
   *  #linprec_->ApplyInverse() since the iterative application is handled by
   *  the Richardson iteration.
   *
   *  \author mayr.mt \date 04/2015
   */
  Teuchos::RCP<Teuchos::ParameterList> CreateIFPACKParameterList() const;

  //! preconditioner
  Teuchos::RCP<LINALG::SOLVER::PreconditionerType> linprec_;

}; // class NlnOperatorRichardson
} // namespace NLNSOL

#endif /* NLN_OPERATOR_RICHARDSON_H */
