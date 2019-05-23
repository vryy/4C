/*----------------------------------------------------------------------------*/
/*!

\brief Wrapper around a linear preconditioner from BACI

\level 3

\maintainer Matthias Mayr
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
  }  // namespace SOLVER
}  // namespace LINALG

namespace NLNSOL
{
  class LineSearchBase;
}  // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnOperatorLinPrec
   *  \brief Wrapper around a linear preconditioner from BACI
   *
   *  \sa NlnOperatorBase
   *
   *  \author mayr.mt \date 10/2014
   */
  class NlnOperatorLinPrec : public NlnOperatorBase
  {
   public:
    //! @name Construction
    //@{

    //! Constructor (empty)
    NlnOperatorLinPrec();

    //! Destructor
    virtual ~NlnOperatorLinPrec() {}

    //@}

    //! @name Setup
    //@{

    //! Setup of the preconditioner object
    virtual void Setup();

    /*! \brief Rebuild preconditioner
     *
     *  No action required.
     *
     *  \author mayr.mt \date 03/2015
     */
    virtual void RebuildPrec(){/* no action required */};

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Apply the preconditioner
     *
     *  We use a linear preconditioner here, that usually solves the system
     *  \f$Ax=b\f$. Though, when reformulating it to an iterative process
     *  \f$A\Delta x = b-Ax_i = -f_i\f$, we can directly use the evaluated
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
    virtual int ApplyInverse(const Epetra_MultiVector& f,  ///< residual in mechanical sense
        Epetra_MultiVector& x                              ///< solution
        ) const;

    //! @name Attribute access functions
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const { return "NLNSOL::NlnOperatorLinPrec"; }

    //! Returns a character string describing the operator briefly
    virtual const char* LabelShort() const { return "LinPrec"; }

    //@}

   protected:
   private:
    //! @name Search Direction
    //@{

    //! Compute line search parameter (=step length)
    virtual void ComputeStepLength(const Epetra_MultiVector& x,  ///< current solution
        const Epetra_MultiVector& f,                             ///< residual at current solution
        const Epetra_MultiVector& inc,                           ///< search direction
        const double fnorm2,                                     ///< current L2-norm of residual
        double& lsparam,  ///< line search parameter (to be determined)
        bool& suffdecr    ///< achieved sufficient decrease? (to be determined)
        ) const;

    //@}

    //! preconditioner
    Teuchos::RCP<LINALG::SOLVER::PreconditionerType> linprec_;

  };  // class NlnOperatorIfpack
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_LIN_PREC_H */
