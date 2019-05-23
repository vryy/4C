/*----------------------------------------------------------------------------*/
/*!

\brief Nonlinear preconditioned Richardson iterations with damping

\level 3

\maintainer Matthias Mayr
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
  /*! \class NlnOperatorRichardson
   *  \brief Richardson iteration with relaxation and preconditioning
   *
   *  This algorithm implements a Richardson iteration with any nonlinear operator
   *  from the NLNSOL framework (denoted by \f$M^{-1}\f$. We allow for a
   *  user-defined relaxation parameter \f$\omega\f$. The iterative procedure
   *  reads
   *  \f[
   *  x_{i+1} = x_i - \omega M^{-1}F(x_i)
   *  \f]
   *
   *  \note In case of IFPACK preconditioners, set the number of sweeps to 1
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
    virtual ~NlnOperatorRichardson() {}

    //@}

    //! @name Setup
    //@{

    //! Setup of the preconditioner object
    virtual void Setup();

    /*! \brief Rebuild preconditioner
     *
     *  \author mayr.mt \date 04/2015
     */
    virtual void RebuildPrec();

    //@}

    //! @name Mathematical functions
    //@{

    /*! \brief Apply the preconditioner
     *
     *  Perform a given number of relaxed and preconditioned Richardson
     *  iterations.
     *
     *  \return Integer error code, set to 0 if successful.
     *
     *  \author mayr.mt \date 04/2015
     */
    virtual int ApplyInverse(const Epetra_MultiVector& f,  ///< residual in mechanical sense
        Epetra_MultiVector& x                              ///< solution
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
    //! @name Setup
    //@{

    //! Setup nonlinear preconditioner
    virtual void SetupPreconditioner();

    //@}

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
    Teuchos::RCP<NLNSOL::NlnOperatorBase> nlnprec_;

  };  // class NlnOperatorRichardson
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_RICHARDSON_H */
