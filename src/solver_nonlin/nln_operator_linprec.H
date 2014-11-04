/*----------------------------------------------------------------------------*/
/*!
 \file nln_operator_linprec.H

 \brief Wrapper around a linear preconditioner from BACI

 <pre>
 Maintainer: Matthias Mayr
             mayr@mhpc.mw.tum.de
             089 - 289-10362
 </pre>
 */

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_LIN_PREC_H
#define NLN_OPERATOR_LIN_PREC_H

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
/*! \brief Wrapper around a linear preconditioner from BACI
 *
 *  \sa NlnOperatorBase
 *
 *  \author mayr.mt \date 10/2014
 */
class NlnOperatorLinPrec: public NlnOperatorBase
{
public:

  //! @name Construction
  //@{

  //! Constructor (empty)
  NlnOperatorLinPrec();

  //! Destructor
  virtual ~NlnOperatorLinPrec()
  {
  }

  //@}

  //! @name Setup
  //@{

  //! Setup of the preconditioner object
  virtual void Setup();

  //@}

  //! @name Mathematical functions
  //@{

  /*! \brief Apply the preconditioner
   *
   *  We use a linear preconditioner here, that usually solves the system
   *  \f$Ax=b\f$. Though, when reformulating it to an interative process
   *  \f$A\Delta x = b-Ax_i = -f_i\f$, we can directly use the evalutated
   *  residual coming from NlnProblem::ComputeF().
   *
   *  To avoid an unnecessary scaling by \f$-1.0\f$, we solve for
   *  \f$A\Delta x = -f_i\f$ and, thus, the update has to read
   *  \f$x_{i+1}=x_i-\Delta x\f$.
   *
   *  \return Integer error code, set to 0 if successful.
   *
   *  \author mayr.mt \date 08/2014
   */
  virtual int ApplyInverse(const Epetra_MultiVector& f, ///< residual (dummy vector to satisfy inheritance from Epetra_Operator)
      Epetra_MultiVector& x ///< solution
      ) const;

  //! @name Attribute access functions
  //@{

  //! Returns a character string describing the operator
  virtual const char* Label() const { return "NLNSOL::NlnOperatorLinPrec"; }

  //@}

protected:

private:

  //! @name Search Direction
  //@{

  //! Compute line search parameter (=step length)
  virtual const double ComputeStepLength(const Epetra_MultiVector& x, ///< current solution
      const Epetra_MultiVector& inc, ///< search direction
      double fnorm2 ///< current L2-norm of residual
      ) const;

  //@}

  //! preconditioner
  Teuchos::RCP<LINALG::SOLVER::PreconditionerType> linprec_;

}; // class NlnOperatorIfpack
} // namespace NLNSOL

#endif /* NLN_OPERATOR_LIN_PREC_H */
